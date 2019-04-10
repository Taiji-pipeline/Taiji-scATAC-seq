{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Quantification
    ( countTagsMerged
    ) where

import Bio.Utils.BitVector
import Bio.Data.Experiment
import qualified Data.ByteString.Char8 as B
import Conduit
import Bio.Data.Bed.Types
import Bio.Data.Bed
import qualified Data.Vector.Unboxed as U
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

-- | Count tags on the merged sample.
countTagsMerged :: SCATACSeqConfig config
                => SCATACSeq S (File '[] 'Bam, a)
                -> WorkflowConfig config (SCATACSeq S (File '[] 'Bed))
countTagsMerged input = do
    genome <- asks (fromJust . _scatacseq_genome_index)
    chrSize <- liftIO $ withGenome genome $ return . getChrSizes
    dir <- asks ((<> "/ReadCount") . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_readcount.bed" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . (\(fl,_) -> do
        runConduit $ countTags fl chrSize 5000 .| writeBed output
        return $ emptyFile & location .~ output )

countTags :: File '[] 'Bam
          -> [(B.ByteString, Int)]
          -> Int  -- ^ resolution
          -> ConduitT i BED IO ()
countTags bam chrSize res = do
    header <- liftIO $ getBamHeader $ bam^.location
    (rc, _) <- liftIO $ runResourceT $ runConduit $
        streamBam (bam^.location) .| bamToBedC header .| 
        countTagsBinBed res beds
    forM_ (zip chrSize rc) $ \((chr, _), vec) -> flip U.imapM_ vec $ \i x -> 
            when (x > 0) $ yield $ BED chr (i * res) ((i+1) * res)
                Nothing (Just $ fromIntegral (x::Int)) Nothing
  where
    beds = map (\(chr, n) -> BED3 chr 0 n) chrSize

    {-
mkCellByBinMatrix :: SCATACSeqConfig config
                => SCATACSeq S (File '[] 'Bam)
                -> WorkflowConfig config (SCATACSeq S (File '[] 'Bed))
                -}

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