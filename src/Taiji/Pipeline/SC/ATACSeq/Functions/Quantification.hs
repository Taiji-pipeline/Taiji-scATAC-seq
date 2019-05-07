{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE GADTs #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Quantification
    ( mkCutSiteIndex
    , getBins
    , mkCountMatrix
    ) where

import Bio.Data.Experiment
import qualified Data.ByteString.Char8 as B
import Data.Conduit.List (groupBy)
import Data.ByteString.Lex.Integral (packDecimal)
import Data.Function (on)
import Data.Conduit.Internal (zipSinks)
import Conduit
import Control.Monad.ST
import Bio.Data.Bed.Types
import Data.Conduit.Zlib (gzip)
import Bio.Data.Bed
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Unboxed.Mutable as UM
import qualified Data.IntervalMap.Strict      as IM
import Control.Monad
import Data.Singletons.Prelude (Elem)
import           Bio.Pipeline.Utils
import Control.Lens
import Control.Monad.Reader (asks)
import Text.Printf (printf)
import qualified Data.Text as T
import Data.Maybe
import Control.DeepSeq (force)
import Bio.Data.Bed.Utils
import Scientific.Workflow
import Bio.Seq.IO (withGenome, getChrSizes)

import Taiji.Pipeline.SC.ATACSeq.Types
import Taiji.Pipeline.SC.ATACSeq.Functions.Utils

-- | Create the cut site map for every cell.
mkCutSiteIndex :: SCATACSeqConfig config
               => SCATACSeq S (File '[NameSorted, Gzip] 'Bed)
               -> WorkflowConfig config (SCATACSeq S (File '[] 'Other))
mkCutSiteIndex input = do
    dir <- asks ((<> "/CutSiteIndex") . _scatacseq_output_dir) >>= getPath
    genome <- fromJust <$> asks _scatacseq_genome_index 
    let output = printf "%s/%s_rep%d.csidx" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . (\fl -> do
        chrs <- withGenome genome $ return . map fst . getChrSizes
        mkCutSiteIndex_ output fl chrs
        return $ emptyFile & location .~ output )
{-# INLINE mkCutSiteIndex #-}

mkCutSiteIndex_ :: FilePath
                -> File '[NameSorted, Gzip] 'Bed
                -> [B.ByteString]
                -> IO ()
mkCutSiteIndex_ output input chrs = createCutSiteIndex output chrs $
    streamBedGzip (input^.location) .| groupBy ((==) `on` (^.name)) .|
    mapC (\x -> (fromJust $ head x ^. name, x))
{-# INLINE mkCutSiteIndex_ #-}

-- | Get bins
getBins :: (Elem 'Gzip tags ~ 'True, SCATACSeqConfig config)
        => SCATACSeq S (File tags 'Bed)
        -> WorkflowConfig config
            (SCATACSeq S (File tags 'Bed, File tags 'Bed, Int))
getBins input = do
    genome <- asks (fromJust . _scatacseq_genome_index)
    chrSize <- liftIO $ withGenome genome $ return . getChrSizes
    dir <- asks ((<> "/ReadCount") . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_bin_idx.bed.gz" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . (\fl -> do
        (n, rc) <- binCount fl chrSize 5000
        let lo = fromIntegral n * 0.001
            hi = fromIntegral n * 0.95
        runResourceT $ runConduit $
            filterBin lo hi 5000 (zip (map fst chrSize) rc) .|
            sinkFileBedGzip output
        return $ (fl, emptyFile & location .~ output, n) )

-- | Divide the genome into bins and count the tags.
binCount :: Elem 'Gzip tags ~ 'True
         => File tags 'Bed    -- ^ Tags
         -> [(B.ByteString, Int)]  -- ^ chromosome sizes
         -> Int  -- ^ resolution
         -> IO (Int, [U.Vector Int])
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

-- | Make the read count matrix.
mkCountMatrix :: (Elem 'Gzip tags ~ 'True, SCATACSeqConfig config)
             => SCATACSeq S (File tags 'Bed, File tags 'Bed, Int)
             -> WorkflowConfig config (SCATACSeq S (File tags 'Other))
mkCountMatrix input = do
    dir <- asks ((<> "/ReadCount") . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_readcount.txt.gz" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . (\(tagFl, regionFl, nCell) -> do
        regions <- runResourceT $ runConduit $
            streamBedGzip (regionFl^.location) .| sinkList :: IO [BED3]
        let bedTree = bedToTree undefined $ zip regions [0::Int ..]
            nBin = length regions
            header = B.pack $ printf "Sparse matrix: %d x %d" nCell nBin
        runResourceT $ runConduit $ streamBedGzip (tagFl^.location) .|
            groupBy ((==) `on` (^.name)) .|
            mapC (tagCountPerCell bedTree nBin) .|
            (yield header >> mapC (uncurry showSparseVector)) .|
            unlinesAsciiC .| gzip .| sinkFile output
        return $ emptyFile & location .~ output )

tagCountPerCell :: BEDTree Int -> Int -> [BED] -> (B.ByteString, U.Vector Int)
tagCountPerCell regions n input = (cellBc, rc)
  where
    rc = U.create $ do
        vec <- UM.replicate n 0
        forM_ input $ \x -> forM_ (IM.elems $ intersecting regions $ toSite x) $
            UM.unsafeModify vec (+1)
        return vec
    cellBc = fromJust $ head input^.name
    toSite bed = case bed^.strand of
        Just False -> BED3 (bed^.chrom) (bed^.chromEnd - 1) (bed^.chromEnd)
        _ -> BED3 (bed^.chrom) (bed^.chromStart) (bed^.chromStart + 1)

showSparseVector :: B.ByteString -> U.Vector Int -> B.ByteString
showSparseVector cellBc = B.intercalate "\t" . (cellBc:) . map f . U.toList .
    U.imapMaybe (\i v -> if v /= 0 then Just (i,v) else Nothing)
  where
    f (i,v) = fromJust (packDecimal i) <> "," <> fromJust (packDecimal v)
{-# INLINE showSparseVector #-}