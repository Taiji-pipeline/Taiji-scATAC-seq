{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Quantification
    ( mkCutSiteIndex
    , countTagsMerged
    ) where

import Bio.Utils.BitVector
import Bio.Data.Experiment
import qualified Data.ByteString.Char8 as B
import Data.Conduit.List (groupBy)
import Data.Function (on)
import Conduit
import Bio.Data.Bed.Types
import Bio.Data.Bed
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Unboxed.Mutable as UM
import qualified Data.IntervalMap.Strict      as IM
import Bio.Data.Bam
import Control.Monad
import           Bio.Pipeline.Utils
import Control.Lens
import Control.Monad.Reader (asks)
import Text.Printf (printf)
import qualified Data.Text as T
import Data.Maybe
import Bio.Data.Bed.Utils
import Scientific.Workflow
import Bio.Seq.IO (withGenome, getChrSizes)

import Taiji.Pipeline.SC.ATACSeq.Types
import Taiji.Pipeline.SC.ATACSeq.Functions.Utils

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

mkCutSiteIndex_ :: FilePath
                -> File '[NameSorted, Gzip] 'Bed
                -> [B.ByteString]
                -> IO ()
mkCutSiteIndex_ output input chrs = createCutSiteIndex output chrs $
    streamBedGzip (input^.location) .| groupBy ((==) `on` (^.name)) .|
    mapC (\x -> (fromJust $ head x ^. name, x))

-- | Count tags on the merged sample.
countTagsMerged :: SCATACSeqConfig config
                => SCATACSeq S (File tags 'Bam, a, Int)
                -> WorkflowConfig config (SCATACSeq S (File tags 'Bed))
countTagsMerged input = do
    genome <- asks (fromJust . _scatacseq_genome_index)
    chrSize <- liftIO $ withGenome genome $ return . getChrSizes
    dir <- asks ((<> "/ReadCount") . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_readcount.bed" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . (\(fl,_,nCells) -> do
        let cutoff = max 5 $ nCells `div` 100
        runResourceT $ runConduit $ countTags fl chrSize 5000 cutoff .|
            sinkFileBed output
        return $ emptyFile & location .~ output )

countTags :: File tags 'Bam
          -> [(B.ByteString, Int)]
          -> Int  -- ^ resolution
          -> Int  -- ^ cutoff
          -> ConduitT i BED (ResourceT IO) ()
countTags bam chrSize res cutoff = do
    header <- liftIO $ getBamHeader $ bam^.location
    (rc, _) <- liftIO $ runResourceT $ runConduit $
        streamBam (bam^.location) .| bamToBedC header .| countTagsBinBed res beds
    forM_ (zip chrSize rc) $ \((chr, _), vec) -> flip U.imapM_ vec $ \i x -> 
            when (x >= cutoff) $ yield $ BED chr (i * res) ((i+1) * res)
                Nothing (Just $ fromIntegral (x::Int)) Nothing
  where
    beds = map (\(chr, n) -> BED3 chr 0 n) chrSize

tagCountPerCell :: Monad m
                => [BED3] -- ^ Open regions
                -> ConduitT [BED] (U.Vector Int) m ()
tagCountPerCell regions = mapC $ \input -> U.create $ do
    vec <- UM.replicate n 0
    forM_ input $ \x -> forM_ (IM.elems $ intersecting bedTree x) $
        UM.unsafeModify vec (+1)
    return vec
  where
    bedTree = bedToTree undefined $ zip regions [0::Int ..]
    n = length regions

{-
mkBaseMap :: FilePath   -- ^ Bam file
          -> [(B.ByteString, Int)]
          -> Int  -- ^ resolution
          -> IO (M.HashMap B.ByteString (U.Vector Bool))
mkBaseMap input chrSize res = runResourceT (
    runConduit $ streamBam input .| bamToBedC .| baseMap chrSize) >>=
    (\BaseMap x -> fmap (downSample res) x) 

-- | Down-sample basemap to given resolution by summing.
downSample :: Int    -- ^ Resolution, must be the multiple of 8.
           -> BitVector
           -> U.Vector Int
downSample res (BitVector _ v) = U.generate n $ \i ->
    U.sum $ U.slice (i * step) step v
  where
    step = res `div` 8
    n = if U.length v `mod` step == 0
        then U.length v `div` step
        else U.length v `div` step + 1
{-# INLINE downSample #-}
-}