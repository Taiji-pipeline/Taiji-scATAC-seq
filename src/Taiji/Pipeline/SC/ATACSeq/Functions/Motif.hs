{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE PartialTypeSignatures #-}

module Taiji.Pipeline.SC.ATACSeq.Functions.Motif
    ( getOpenRegion
    , findMotifsPre
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
import Data.Binary
import Control.Monad.Reader (asks)
import           Bio.Motif                     hiding (score)
import Data.Maybe
import Conduit
import Control.Monad
import qualified Data.HashMap.Strict as M
import qualified Data.ByteString.Char8 as B
import qualified Bio.Utils.BitVector as BV
import           System.IO
import           Bio.Seq.IO
import           Data.Default (def)
import           Shelly                        (fromText, mkdir_p, shelly,
                                                test_f)
import           System.FilePath               (takeDirectory)

import Taiji.Pipeline.SC.ATACSeq.Types

findMotifsPre :: SCATACSeqConfig config
              => Double
              -> File '[Gzip] 'Bed
              -> WorkflowConfig config
                  [(B.ByteString, File '[Gzip] 'Bed, File '[Gzip] 'Other)]
findMotifsPre p region = do
    motifFile <- fromMaybe (error "Motif file is not specified!") <$>
        asks _scatacseq_motif_file
    genome <- asks (fromJust . _scatacseq_genome_index)
    chrs <- liftIO $ withGenome genome $ return . map fst . getChrSizes
    dir <- asks _scatacseq_output_dir >>= getPath . (<> (asDir "/temp"))
    let output = dir ++ "/motif.bin"
    liftIO $ map (mkCutoffMotif def p) <$> readMEME motifFile >>= encodeFile output
    return $ zip3 chrs (repeat region) $ repeat $ location .~ output $ emptyFile

-- | Identify motif binding sites.
findMotifs :: SCATACSeqConfig config
           => (B.ByteString, File '[Gzip] 'Bed, File '[Gzip] 'Other)
           -> WorkflowConfig config (File '[Gzip] 'Bed)
findMotifs (chr, openChromatin, motifFl) = do
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
    dir <- asks _scatacseq_output_dir >>= getPath . (<> (asDir "/TFBS"))
    let output = dir ++ "/" ++ B.unpack chr ++ ".bed.gz"
    liftIO $ withGenome seqIndex $ \g -> do
        motifs <- decodeFile (motifFl^.location) :: IO [CutoffMotif]
        runResourceT $ runConduit $
            (streamBedGzip (openChromatin^.location) :: _ _ BED3 _ _) .|
            filterC (\x -> x^.chrom == chr) .| scanMotif g motifs .|
            sinkFileBedGzip output
        return $ location .~ output $ emptyFile
{-# INLINE findMotifs #-}

-- mkBigBed :: File '[Gzip] 'Bed

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
    runResourceT $ runConduit $ source .|
        (baseMap chrSize >>= baseMapToRegion) .| sinkFileBedGzip output
    return $ emptyFile & location .~ output
{-# INLINE getOpenRegion #-}

baseMapToRegion :: Monad m => BaseMap -> ConduitT i BED3 m ()
baseMapToRegion (BaseMap bm) = forM_ (M.toList bm) $ \(chr, bv) ->
    mkBedWindow 50 chr bv
{-# INLINE baseMapToRegion #-}

-- | Create bed from basemap.
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
{-# INLINE mkBedWindow #-}
