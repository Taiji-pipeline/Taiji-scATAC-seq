{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE PartialTypeSignatures #-}

module Taiji.Pipeline.SC.ATACSeq.Functions.Motif
    ( getOpenRegion
    , prepDataSet
    , findMotifs
    ) where

import Bio.Data.Experiment
import Control.Lens
import qualified Data.Text as T
import Scientific.Workflow
import Bio.Seq.IO (withGenome, getChrSizes)
import Bio.Data.Bed.Utils
import           Bio.Pipeline
import           Bio.Pipeline.Instances ()
import Bio.Data.Bed
import Data.List.Split (chunksOf)
import Control.Monad.Reader (asks)
import           Bio.Motif                     hiding (score)
import Data.Maybe
import Conduit
import Control.Monad
import qualified Data.HashMap.Strict as M
import qualified Data.ByteString.Char8 as B
import qualified Bio.Utils.BitVector as BV
import           System.IO
import           System.IO.Temp                (emptyTempFile)
import           Bio.Seq.IO
import           Data.Default (def)
import           Shelly                        (fromText, mkdir_p, shelly,
                                                test_f)
import           System.FilePath               (takeDirectory)

import Taiji.Pipeline.SC.ATACSeq.Types

getOpenRegion :: SCATACSeqConfig config
              => [SCATACSeq S (File '[NameSorted, Gzip] 'Bed)]
              -> WorkflowConfig config (File '[Gzip] 'Bed)
getOpenRegion inputs = do
    dir <- asks _scatacseq_output_dir >>= getPath
    genome <- asks (fromJust . _scatacseq_genome_index)
    chrSize <- liftIO $ withGenome genome $ return . getChrSizes
    let output = dir <> "/accessible_region.bed.gz"
        source = forM_ inputs $ \input ->
            streamBedGzip $ input^.replicates._2.files.location
    runResourceT $ runConduit $ source .| getOpenRegion_ chrSize .|
        sinkFileBedGzip output
    return $ emptyFile & location .~ output

getOpenRegion_ :: PrimMonad m
               => [(B.ByteString, Int)]  -- ^ Chromosomes and their sizes
               -> ConduitT BED BED3 m ()
getOpenRegion_ chrs = baseMap chrs >>= baseMapToRegion

baseMapToRegion :: Monad m => BaseMap -> ConduitT i BED3 m ()
baseMapToRegion (BaseMap bm) = forM_ (M.toList bm) $ \(chr, bv) ->
    mkBedWindow 50 chr bv

mkBedWindow :: Monad m 
            => Int   -- ^ half window size
            -> B.ByteString
            -> BV.BitVector -> ConduitT i BED3 m ()
mkBedWindow window chr bv = go (0,0) 0
  where
    go (lo, hi) i
        | i >= n = yield $ asBed chr lo hi
        | not (bv BV.! i) = go (lo, hi) $! i + 1
        | cur_lo <= hi = go (lo, cur_hi) $! i + 1
        | otherwise = yield (asBed chr lo hi) >> go (cur_lo, cur_hi) (i + 1)
      where
        (cur_lo, cur_hi) = getRegion i
    n = BV.size bv
    getRegion x = (max 0 $ x - window, min n $ x + 1 + window)

prepDataSet :: SCATACSeqConfig config
            => File '[Gzip] 'Bed
            -> WorkflowConfig config
                (ContextData (File '[Gzip] 'Bed) [[Motif]])
prepDataSet region = do
    motifFile <- fromMaybe (error "Motif file is not specified!") <$>
        asks _scatacseq_motif_file
    motifs <- liftIO $ readMEME motifFile
    return $ ContextData region $ chunksOf 100 motifs

findMotifs :: SCATACSeqConfig config
           => Double     -- ^ p value
           -> ContextData (File '[Gzip] 'Bed) [Motif]
           -> WorkflowConfig config (File '[Gzip] 'Bed)
findMotifs p (ContextData openChromatin motifs) = do
    -- Generate sequence index
    genome <- asks ( fromMaybe (error "Genome fasta file was not specified!") .
        _scatacseq_genome_fasta )
    seqIndex <- asks ( fromMaybe (error "Genome index file was not specified!") .
        _scatacseq_genome_index )
    fileExist <- liftIO $ shelly $ test_f $ fromText $ T.pack seqIndex
    liftIO $ if fileExist
        then hPutStrLn stderr "Sequence index exists. Skipped."
        else do
            shelly $ mkdir_p $ fromText $ T.pack $ takeDirectory seqIndex
            hPutStrLn stderr "Generating sequence index"
            mkIndex [genome] seqIndex
    dir <- asks _scatacseq_output_dir >>= getPath . (<> (asDir "/TFBS/"))
    liftIO $ withGenome seqIndex $ \g -> do
        output <- emptyTempFile dir "motif_sites_part.bed.gz"
        runResourceT $ runConduit $
            (streamBedGzip (openChromatin^.location) :: _ _ BED3 _ _) .|
            motifScan g motifs def p .| getMotifScore g motifs def .|
            getMotifPValue (Just (1 - p * 10)) motifs def .| sinkFileBedGzip output
        return $ location .~ output $ emptyFile