{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE GADTs #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Feature.Window
    ( getWindows
    , mkWindowMat
    ) where

import qualified Data.ByteString.Char8 as B
import Data.Conduit.List (groupBy)
import Data.Conduit.Internal (zipSinks)
import Control.Monad.ST
import Bio.Data.Bed.Types
import Bio.Data.Bed
import qualified Data.Vector.Unboxed as U
import Data.Singletons.Prelude (Elem)
import           Bio.Pipeline.Utils
import qualified Data.Text as T
import Control.DeepSeq (force)
import Bio.Data.Bed.Utils
import Bio.Seq.IO (withGenome, getChrSizes)

import Taiji.Prelude hiding (groupBy)
import Taiji.Utils
import Taiji.Pipeline.SC.ATACSeq.Types
import Taiji.Pipeline.SC.ATACSeq.Functions.Utils

-- | Get candidate bins.
getWindows :: (Elem 'Gzip tags ~ 'True, SCATACSeqConfig config)
           => FilePath
           -> SCATACSeq S (File tags 'Bed)
           -> ReaderT config IO
              (SCATACSeq S (File tags 'Bed, File tags 'Bed, Int))
getWindows prefix input = do
    res <- fromIntegral <$> asks _scatacseq_window_size
    genome <- asks (fromJust . _scatacseq_genome_index)
    chrSize <- liftIO $ withGenome genome $ return . getChrSizes
    dir <- asks ((<> asDir prefix) . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_window_idx.bed.gz" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ ( \fl -> do
        (n, rc) <- liftIO $ binCount fl chrSize res
        let windows = forM_ (zip (map fst chrSize) rc) $ \(chr, vec) ->
                flip U.imapM_ vec $ \i x -> when (x > 0) $ yield $
                    BED chr (i * res) ((i+1) * res) Nothing (Just x) Nothing
        asks _scatacseq_blacklist >>= \case
            Nothing -> liftIO $ do
                runResourceT $ runConduit $ windows .| sinkFileBedGzip output
            Just blacklist -> liftIO $ do
                blackRegions <- readBed blacklist :: IO [BED3]
                let bedTree = bedToTree const $ map (\x -> (x, ())) blackRegions
                runResourceT $ runConduit $ windows .|
                    filterC (not . isIntersected bedTree) .| sinkFileBedGzip output
        return $ (fl, emptyFile & location .~ output, n) )

-- | Make the read count matrix.
mkWindowMat :: (Elem 'Gzip tags ~ 'True, SCATACSeqConfig config)
               => FilePath
               -> SCATACSeq S (File tags 'Bed, File tags 'Bed, Int)
               -> ReaderT config IO (SCATACSeq S (File tags 'Matrix))
mkWindowMat dir input = do
    let output = printf "%s/%s_rep%d_window.mat.gz" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . (\(tagFl, regionFl, nCell) -> do
        regions <- runResourceT $ runConduit $
            streamBedGzip (regionFl^.location) .| sinkList :: IO [BED3]
        runResourceT $ runConduit $ streamBedGzip (tagFl^.location) .|
            groupCells .| mkCountMat regions .|
            sinkRows nCell (length regions) (fromJust . packDecimal) output
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