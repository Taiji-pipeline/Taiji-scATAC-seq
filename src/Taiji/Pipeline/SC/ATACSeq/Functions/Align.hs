{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE RecordWildCards #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Align
    ( tagAlign
    , mkIndices
    , filterBamSort
    , deDuplicates
    , filterCell
    , mkBigWig
    ) where

import Bio.Data.Bed
import Bio.Data.Bed.Types (BED(..), BED3(..))
import Bio.Data.Bam
import           Bio.HTS
import Data.Conduit.Internal (zipSinks)
import Control.Monad.State.Strict
import Data.Conduit.List (groupBy)
import qualified Data.Text as T
import Shelly hiding (FilePath)
import qualified Data.HashSet as S
import qualified Data.ByteString.Char8 as B
import System.IO
import           Bio.Seq.IO (mkIndex, withGenome, getChrSizes)
import           System.FilePath               (takeDirectory)

import Taiji.Prelude hiding (groupBy)
import Taiji.Pipeline.SC.ATACSeq.Types
import Taiji.Pipeline.SC.ATACSeq.Functions.Utils
import Taiji.Pipeline.SC.ATACSeq.Functions.QC

mkIndices :: SCATACSeqConfig config => [a] -> ReaderT config IO ()
mkIndices input
    | null input = return ()
    | otherwise = do
        genome <- asks $ fromJust . _scatacseq_genome_fasta

        -- Generate BWA index
        dir <- asks (fromJust . _scatacseq_bwa_index)
        _ <- liftIO $ bwaMkIndex genome dir

        -- Generate genome index
        seqIndex <- asks $ fromMaybe
            (error "Genome index file was not specified!") .
            _scatacseq_genome_index
        fileExist <- liftIO $ shelly $ test_f $ fromText $ T.pack seqIndex
        liftIO $ if fileExist
            then hPutStrLn stderr "Sequence index exists. Skipped."
            else do
                shelly $ mkdir_p $ fromText $ T.pack $ takeDirectory seqIndex
                hPutStrLn stderr "Generating sequence index"
                mkIndex [genome] seqIndex
        return ()

tagAlign :: SCATACSeqConfig config
         => SCATACSeq S ( Either (SomeTags 'Fastq)
                (SomeTags 'Fastq, SomeTags 'Fastq) )
         -> ReaderT config IO ( SCATACSeq S
                (Either (File '[] 'Bam) (File '[PairedEnd] 'Bam)) )
tagAlign input = do
    dir <- asks ((<> "/Bam") . _scatacseq_output_dir) >>= getPath
    idx <- asks (fromJust . _scatacseq_bwa_index)
    let output = printf "%s/%s_rep%d.bam" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . ( \fl -> case fl of
        Left f ->
            let f' = fromSomeTags f :: File '[] 'Fastq
            in bwaAlign output idx (Left f') $ defaultBWAOpts & bwaCores .~ 4
        Right (f1,f2) ->
            let f1' = fromSomeTags f1 :: File '[] 'Fastq
                f2' = fromSomeTags f2 :: File '[] 'Fastq
            in bwaAlign output idx (Right (f1', f2')) $
                defaultBWAOpts & bwaCores .~ 4
        )

filterBamSort :: SCATACSeqConfig config
              => SCATACSeq S (Either (File '[] 'Bam) (File '[PairedEnd] 'Bam))
              -> ReaderT config IO ( SCATACSeq S ( Either
                  (File '[NameSorted] 'Bam)
                  (File '[NameSorted, PairedEnd] 'Bam) ))
filterBamSort input = do
    dir <- asks ((<> "/Bam") . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_srt.bam" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . either
        (fmap Left . fun output) (fmap Right . filterBam "./" output)
  where
    fun output x = withTemp (Just "./") $ \f ->
        filterBam "./" f x >>= sortBamByName "./" output

deDuplicates :: SCATACSeqConfig config
             => SCATACSeq S ( Either
                 (File '[NameSorted] 'Bam)
                 (File '[NameSorted, PairedEnd] 'Bam) )
             -> ReaderT config IO
                 (SCATACSeq S ( File '[NameSorted, Filtered] 'Bam  -- ^ filtered reads
                              , File '[NameSorted] 'Bam        -- ^ mito reads
                              , File '[] 'Tsv ))               -- ^ Stat
deDuplicates input = do
    dir <- asks _scatacseq_output_dir >>= getPath . (<> (asDir "/Bam"))
    tss <- asks _scatacseq_annotation >>= liftIO . readTSS . fromJust
    let outputBam = printf "%s/%s_rep%d_srt_filt.bam" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
        outputStat = printf "%s/%s_rep%d_stat.txt" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
        outputMito = printf "%s/%s_rep%d_mito.bam" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
        f fl = do
            header <- getBamHeader fl
            _ <- runResourceT $ runConduit $ streamBam fl .|
                groupBy ((==) `on` (extractBarcode . queryName)) .|
                mapC (filterReads header tss) .| zipSinks
                    ( zipSinks
                        (concatMapC (fst . fst) .| sinkBam outputBam header)
                        (concatMapC (snd . fst) .| sinkBam outputMito header) )
                    (mapC (showStat . snd) .| unlinesAsciiC .| sinkFile outputStat)
            return ( location .~ outputBam $ emptyFile
                   , location .~ outputMito $ emptyFile
                   , location .~ outputStat $ emptyFile )
    input & replicates.traverse.files %%~ liftIO . f . either
        (^.location) (^.location)

-- | Remove duplicated and chrM reads. 
filterReads :: BAMHeader
            -> BEDTree (Int, Bool)   -- ^ TSS
            -> [BAM]
            -> (([BAM], [BAM]), Stat)
filterReads hdr tss bam = runState filterFn $ Stat bc 0 0 0 0 0
  where
    filterFn = do
        (tags, mito) <- rmAbnoramlFragment bam >>= rmDup >>= rmChrM hdr
        tssEnrichment tss hdr tags
        totalReads tags
        return (tags, mito)
    bc = extractBarcode $ queryName $ head bam
{-# INLINE filterReads #-}

-- | Filter QC-failed cells and convert bam to bed with TN5 correction.
-- The BED interval of the fragment is obtained by adjusting the BAM alignment
-- interval of the sequenced read-pair. The start of the interval is moved
-- forward by 4bp from a left-most alignment position and backward 5bp from the
-- right-most alignment position. The transposase cuts the two DNA strands with
-- a 9bp overhang, and adjusted positions represent the center point between
-- these cuts; this position is recorded as a cut site that represents a
-- chromatin accessibility event.
filterCell :: SCATACSeqConfig config
           => (SCATACSeq S ( File '[NameSorted, Filtered] 'Bam 
                           , File '[NameSorted] 'Bam        
                           , File '[] 'Tsv ))               
           -> ReaderT config IO (SCATACSeq S (File '[NameSorted, Gzip] 'Bed))
filterCell input = do
    dir <- asks _scatacseq_output_dir >>= getPath . (<> (asDir "/Bed"))
    let output = printf "%s/%s_rep%d_srt_filt.bed.gz" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . (\(bamFl, _, statFl) -> do
        stats <- readStats $ statFl^.location
        let cells = S.fromList $ map _barcode $ filter (passedQC 1) stats
        header <- getBamHeader $ bamFl^.location
        runResourceT $ runConduit $ streamBam (bamFl^.location) .|
            mapC (toBed header) .| groupBy ((==) `on` (^.name)) .|
            filterC ((`S.member` cells) . fromJust . (^.name) . head) .|
            concatC .| sinkFileBedGzip output
        return $ location .~ output $ emptyFile )
  where
    toBed header bam
        | tLen bam > 0 = BED chr (start+4) end nm sc str
        | otherwise = BED chr start (end-5) nm sc str
      where
        chr = fromJust $ refName header bam 
        start = startLoc bam
        end = endLoc bam
        nm = Just $ extractBarcode $ queryName bam
        str = Just $ not $ isRev bam
        sc = Just $ fromIntegral $ mapq bam
{-# INLINE filterCell #-}

mkBigWig :: SCATACSeqConfig config
         => FilePath   -- output
         -> File tags 'Bed
         -> ReaderT config IO ()
mkBigWig output input = do
    seqIndex <- getGenomeIndex
    liftIO $ withTempDir (Just "./") $ \dir -> do
        let tmp1 = dir ++ "/tmp1"
            tmp2 = dir ++ "/tmp2"
            tmpChr = dir ++ "/chr"
        extendBed tmp1 $ input^.location
        shelly $ escaping False $ run_ "sort"
            ["-k", "1,1", "-k2,2n", T.pack tmp1, ">", T.pack tmp2]
        mkBedGraph tmp1 tmp2
        chrSize <- withGenome seqIndex $ return . getChrSizes
        B.writeFile tmpChr $ B.unlines $
            map (\(a,b) -> a <> "\t" <> B.pack (show b)) chrSize
        shelly $ run_ "bedGraphToBigWig" [T.pack tmp1, T.pack tmpChr, T.pack output]
  where
    extendBed out fl = runResourceT $ runConduit $
        streamBedGzip fl .| mapC f .| sinkFileBed out
      where
        f :: BED -> BED3
        f bed = case bed^.strand of
            Just False -> BED3 (bed^.chrom) (bed^.chromEnd - 100) (bed^.chromEnd)
            _ -> BED3 (bed^.chrom) (bed^.chromStart) (bed^.chromStart + 100)

mkBedGraph :: FilePath  -- ^ Output
           -> FilePath  -- ^ Coordinate sorted bed files
           -> IO ()
mkBedGraph output input = runResourceT $ runConduit $ streamBedGzip input .|
    mergeSortedBedWith (splitOverlapped length :: [BED3] -> [(BED3, Int)]) .|
    concatC .| mapC f .| unlinesAsciiC .| sinkFile output
  where
    f (bed, x) = toLine bed <> "\t" <> B.pack (show x)
{-# INLINE mkBedGraph #-}