{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE OverloadedLists #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE GADTs #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Align
    ( tagAlign
    , mkIndices
    , filterNameSortBam
    , sortBedFile
    , deDuplicates
    , getQCMetric
    , filterCell
    ) where

import Bio.Data.Bed
import Bio.Data.Bam
import           Bio.HTS
import Control.Monad.State.Strict
import Data.Conduit.List (groupBy)
import qualified Data.Text as T
import qualified Data.HashSet as S
import           Data.Coerce             (coerce)
import Control.DeepSeq (force)
import Control.Exception (evaluate)
import Shelly (shelly, bash_, bashPipeFail, silently, escaping)

import Taiji.Prelude hiding (groupBy)
import Taiji.Pipeline.SC.ATACSeq.Types
import Taiji.Pipeline.SC.ATACSeq.Functions.Utils
import Taiji.Pipeline.SC.ATACSeq.Functions.QC

mkIndices :: SCATACSeqConfig config
          => [SCATACSeq N [Either SomeFile (SomeFile, SomeFile)]]
          -> ReaderT config IO
              [SCATACSeq N [Either SomeFile (SomeFile, SomeFile)]]
mkIndices [] = return []
mkIndices input = do
    mkGenomeIndex 
    let fq = filter (isFq . either id fst) $ concat $
            input^..folded.replicates.folded.files
    unless (null fq) $ do
        genome <- getGenomeFasta
        -- Generate BWA index
        idx <- asks (fromJust . _scatacseq_bwa_index)
        liftIO (bwaMkIndex genome idx) >> return ()
    return input
  where
    isFq x = getFileType x == Fastq

tagAlign :: SCATACSeqConfig config
         => SCATACSeq S (File '[Demultiplexed, Gzip] 'Fastq, File '[Demultiplexed, Gzip] 'Fastq)
         -> ReaderT config IO ( SCATACSeq S
                (Either (File '[] 'Bam) (File '[PairedEnd] 'Bam)) )
tagAlign input = do
    dir <- asks ((<> "/Bam") . _scatacseq_output_dir) >>= getPath
    idx <- asks (fromJust . _scatacseq_bwa_index)
    let output = printf "%s/%s_rep%d.bam" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . ( \(f1, f2) ->
        let f1' = location .~ (f1^.location) $ emptyFile :: File '[] 'Fastq
            f2' = location .~ (f2^.location) $ emptyFile :: File '[] 'Fastq
        in bwaAlign output idx (Right (f1', f2')) $ defaultBWAOpts & bwaCores .~ 8
        )

filterNameSortBam :: SCATACSeqConfig config
                  => SCATACSeq S (Either (File '[] 'Bam) (File '[PairedEnd] 'Bam))
                  -> ReaderT config IO ( SCATACSeq S ( Either
                      (File '[NameSorted] 'Bam)
                      (File '[NameSorted, PairedEnd] 'Bam) ))
filterNameSortBam input = do
    dir <- asks ((<> "/Bam") . _scatacseq_output_dir) >>= getPath
    tmp <- fromMaybe "./" <$> asks _scatacseq_tmp_dir
    let output = printf "%s/%s_rep%d_srt.bam" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
        fun x = withTemp (Just tmp) $ \f ->
            filterBam tmp f x >>= sortBamByName tmp output
    input & replicates.traverse.files %%~ liftIO . either
        (fmap Left . fun) (fmap Right . filterBam tmp output)

sortBedFile :: ( SCATACSeqConfig config
               , unpair ~ '[NameSorted, Gzip]
               , paired ~ '[NameSorted, PairedEnd, Gzip] )
            => SCATACSeq S (SomeTags 'Bed)
            -> ReaderT config IO (
                 SCATACSeq S (Either (File unpair 'Bed) (File paired 'Bed)) )
sortBedFile input = do
    dir <- asks _scatacseq_output_dir >>= getPath . (<> (asDir "/Bed"))
    tmpdir <- fromMaybe "./" <$> asks _scatacseq_tmp_dir
    let output = printf "%s/%s_rep%d_nsrt_fragments.bed.gz" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
        f (SomeTags fl) = do
            res <- if fl `hasTag` NameSorted
                then if fl `hasTag` Gzip
                    then return $ fl^.location
                    else do
                        shelly $ escaping False $ silently $ bashPipeFail bash_ "gzip"
                            ["-c", ">", T.pack output]
                        return output
                else if fl `hasTag` Gzip
                    then do
                        shelly $ escaping False $ silently $ bashPipeFail bash_ "gzip"
                            [ "-d", "-c", T.pack $ fl^.location, "|"
                            , "sort", "-T", T.pack tmpdir, "-k4,4", "|"
                            , "gzip", "-c", ">", T.pack output ]
                        return output
                    else do
                        shelly $ escaping False $ silently $ bashPipeFail bash_ "sort"
                            [ "-T", T.pack tmpdir, "-k4,4", T.pack $ fl^.location, "|"
                            , "gzip", "-c", ">", T.pack output ]
                        return output
            return $ if fl `hasTag` PairedEnd
                then Right $ location .~ res $ emptyFile
                else Left $ location .~ res $ emptyFile
    input & replicates.traverse.files %%~ liftIO . f

-- | Remove duplicates and convert bam file to fragments.
-- The BED interval of the fragment is obtained by adjusting the BAM alignment
-- interval of the sequenced read-pair. The start of the interval is moved
-- forward by 4bp from a left-most alignment position and backward 5bp from the
-- right-most alignment position. The transposase cuts the two DNA strands with
-- a 9bp overhang, and adjusted positions represent the center point between
-- these cuts; this position is recorded as a cut site that represents a
-- chromatin accessibility event.
deDuplicates :: ( SCATACSeqConfig config
                , unpair ~ '[NameSorted]
                , paired ~ '[NameSorted, PairedEnd]
                , unpair' ~ Insert' Gzip unpair
                , paired' ~ Insert' Gzip paired )
             => SCATACSeq S (Either (File unpair 'Bam) (File paired 'Bam))
             -> ReaderT config IO (
                 SCATACSeq S (Either (File unpair' 'Bed) (File paired' 'Bed)) )
deDuplicates input = do
    dir <- asks _scatacseq_output_dir >>= getPath . (<> (asDir "/Bed"))
    let output = printf "%s/%s_rep%d_nsrt_fragments.bed.gz" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . either
        (fmap (Left . coerce) . f output . SomeFile)
        (fmap Right . f output . SomeFile)
  where
    f output (SomeFile fl) = do
        let pair = fl `hasTag` PairedEnd
        header <- getBamHeader $ fl^.location
        runResourceT $ runConduit $ streamBam (fl^.location) .|
            groupBy ((==) `on` (extractBarcode . queryName)) .|
            concatMapC (\bam -> rmDup pair header (if pair then rmAbnoramlFragment bam else bam)) .|
            sinkFileBedGzip output
        return $ location .~ output $ emptyFile

getQCMetric :: ( SCATACSeqConfig config
               , unpair ~ '[NameSorted, Gzip]
               , paired ~ '[NameSorted, PairedEnd, Gzip] )
             => SCATACSeq S (Either (File unpair 'Bed) (File paired 'Bed))
             -> ReaderT config IO (SCATACSeq S 
                 ( Either (File unpair 'Bed) (File paired 'Bed)
                 , File '[] 'Tsv ))
getQCMetric input = do
    dir <- qcDir
    tss <- getAnnotation >>= liftIO . readTSS
    let output = printf "%s/%s_rep%d_qc.tsv" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . \fl -> do
        q <- either (qc output tss . SomeFile) (qc output tss . SomeFile) fl
        return (fl, q)
  where
    qc output tss (SomeFile fl) = do
        stats <- runResourceT $ runConduit $ streamBedGzip (fl^.location) .|
            groupBy ((==) `on` (^.name)) .| mapMC f .| sinkList
        writeStats output stats
        return $ location .~ output $ emptyFile
      where
        f beds = liftIO $ evaluate $ force $ Stat bc dupRate (Just chrMRate) te num Nothing
          where
            bc = fromJust $ head beds ^. name
            dupRate = case totalReads of
                Nothing -> Nothing
                Just s -> Just $ 1 - fromIntegral num / fromIntegral s 
            (beds', _, chrMRate) = rmChrM beds
            te = tssEnrichment tss beds'
            num = length beds'
            totalReads = fmap (foldl1' (+)) $ sequence $ map (^.score) beds

-- | Filter QC-failed cells and convert fragments into TN5 insertions.
filterCell :: ( SCATACSeqConfig config
              , unpair ~ '[NameSorted, Gzip]
              , paired ~ '[NameSorted, PairedEnd, Gzip] )
           => SCATACSeq S ( Either (File unpair 'Bed) (File paired 'Bed)
                          , File '[] 'Tsv )
           -> ReaderT config IO (SCATACSeq S (File '[NameSorted, Gzip] 'Bed))
filterCell input = do
    dir <- asks _scatacseq_output_dir >>= getPath . (<> (asDir "/Bed"))
    passedQC <- getQCFunction
    let output = printf "%s/%s_rep%d_TN5_insertion.bed.gz" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . ( \(bedFl, qcFl) -> do
        stats <- readStats $ qcFl^.location
        let cells = S.fromList $ map _barcode $ filter passedQC stats
        either (f output cells . SomeFile) (f output cells . SomeFile) bedFl
        return $ location .~ output $ emptyFile )
  where
    f output cells (SomeFile fl) = runResourceT $ runConduit $
        streamBedGzip (fl^.location) .| groupBy ((==) `on` (^.name)) .|
        filterC ((`S.member` cells) . fromJust . (^.name) . head) .|
        concatMapC (concatMap getCutSite . filter notChrM) .| sinkFileBedGzip output
      where
        getCutSite x | fl `hasTag` PairedEnd = [left, right]
                     | otherwise = [x]
          where
            left = let i = x^.chromStart in BED (x^.chrom) i (i+1) (x^.name) Nothing (Just True)
            right = let i = x^.chromEnd - 1 in BED (x^.chrom) i (i+1) (x^.name) Nothing (Just False)
    notChrM x = let chr = x^.chrom in chr /= "chrM" && chr /= "M"
