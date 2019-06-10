{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE GADTs #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Feature.Window
    ( getBins
    , mkWindowMat
    ) where

import qualified Data.ByteString.Char8 as B
import Data.Conduit.List (groupBy)
import qualified Data.HashMap.Strict                  as M
import qualified Data.HashSet as S
import Data.ByteString.Lex.Integral (packDecimal)
import Data.Double.Conversion.ByteString (toShortest)
import           Bio.Utils.Misc (readInt, readDouble)
import Data.Conduit.Internal (zipSinks)
import Control.Arrow (first)
import Control.Monad.ST
import Bio.Data.Bed.Types
import Data.Conduit.Zlib (gzip)
import Bio.Data.Bed
import Bio.RealWorld.GENCODE (readGenes, Gene(..))
import           Data.CaseInsensitive  (mk, original, CI)
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Unboxed.Mutable as UM
import qualified Data.IntervalMap.Strict      as IM
import Data.Singletons.Prelude (Elem)
import           Bio.Pipeline.Utils
import Control.Arrow (second)
import           Data.List.Ordered                    (nubSort)
import qualified Data.Text as T
import Control.DeepSeq (force)
import Bio.Data.Bed.Utils
import Bio.Seq.IO (withGenome, getChrSizes)

import Taiji.Prelude hiding (groupBy)
import Taiji.Pipeline.SC.ATACSeq.Types
import Taiji.Pipeline.SC.ATACSeq.Functions.Utils

-- | Get candidate bins.
getBins :: (Elem 'Gzip tags ~ 'True, SCATACSeqConfig config)
        => SCATACSeq S (File tags 'Bed)
        -> ReaderT config IO
            (SCATACSeq S (File tags 'Bed, File tags 'Bed, Int))
getBins input = do
    genome <- asks (fromJust . _scatacseq_genome_index)
    chrSize <- liftIO $ withGenome genome $ return . getChrSizes
    dir <- asks ((<> "/ReadCount") . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_bin_idx.bed.gz" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . (\fl -> do
        (n, rc) <- binCount fl chrSize res
        let lo = fromIntegral n * 0.001
            hi = fromIntegral n * 1
        runResourceT $ runConduit $
            filterBin lo hi res (zip (map fst chrSize) rc) .|
            sinkFileBedGzip output
        return $ (fl, emptyFile & location .~ output, n) )
  where
    res = 5000

-- | Make the read count matrix.
mkWindowMat :: (Elem 'Gzip tags ~ 'True, SCATACSeqConfig config)
               => SCATACSeq S (File tags 'Bed, File tags 'Bed, Int)
               -> ReaderT config IO (SCATACSeq S (File tags 'Other))
mkWindowMat input = do
    dir <- asks ((<> "/ReadCount") . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_window.txt.gz" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . (\(tagFl, regionFl, nCell) -> do
        regions <- runResourceT $ runConduit $
            streamBedGzip (regionFl^.location) .| sinkList :: IO [BED3]
        runResourceT $ runConduit $ streamBedGzip (tagFl^.location) .|
            groupCells .| mkFeatMat nCell regions .| sinkFile output
        return $ emptyFile & location .~ output )

-- | Divide the genome into bins and count the tags.
binCount :: Elem 'Gzip tags ~ 'True
         => File tags 'Bed    -- ^ Tags
         -> [(B.ByteString, Int)]  -- ^ chromosome sizes
         -> Int  -- ^ resolution
         -> IO (Int, [U.Vector Int])   -- ^ Cell number and count for each chromosome
binCount input chrSize res = runResourceT $ runConduit $
    inputWithGroupCell .| zipSinks lengthC (mapC countTags .| fold)
  where
    inputWithGroupCell = streamBedGzip (input^.location) .|
        groupBy ((==) `on` (^.name))
    countTags beds = map (U.map (\x -> if x > 0 then 1 else x)) $
        fst $ runST $ runConduit $ yieldMany beds .| countTagsBinBed res chrs
    chrs = map (\(chr, n) -> BED3 chr 0 n) chrSize
    fold = await >>= \case
        Nothing -> return []
        Just v -> foldlC f v
    f x y = force $ zipWith (U.zipWith (+)) x y
{-# INLINE binCount #-}

filterBin :: Monad m
          => Double
          -> Double
          -> Int
          -> [(B.ByteString, U.Vector Int)] -> ConduitT i BED m ()
filterBin lo hi res rc = forM_ rc $ \(chr, vec) -> flip U.imapM_ vec $ \i x -> 
    when (fromIntegral x > lo && fromIntegral x < hi) $
        yield $ BED chr (i * res) ((i+1) * res) Nothing (Just x) Nothing
{-# INLINE filterBin #-}