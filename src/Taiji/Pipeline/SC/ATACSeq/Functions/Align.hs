{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE RecordWildCards #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Align
    ( tagAlign
    , mkIndex
    , filterBamSort
    , qualityControl
    ) where

import Bio.Data.Bed
import Bio.Data.Bam
import           Bio.HTS
import Data.Conduit.Internal (zipSinks)
import Control.Monad.State.Strict
import Data.Conduit.List (groupBy)
import qualified Data.Text as T

import Taiji.Prelude hiding (groupBy, frip)
import Taiji.Pipeline.SC.ATACSeq.Types
import Taiji.Pipeline.SC.ATACSeq.Functions.Utils
import Taiji.Pipeline.SC.ATACSeq.Functions.QC

mkIndex :: SCATACSeqConfig config
        => [a] -> ReaderT config IO ()
mkIndex input
    | null input = return ()
    | otherwise = do
        genome <- asks (fromJust . _scatacseq_genome_fasta)
        -- Generate BWA index
        dir <- asks (fromJust . _scatacseq_bwa_index)
        _ <- liftIO $ bwaMkIndex genome dir
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
            in bwaAlign output idx (Left f') $ defaultBWAOpts & bwaCores .~ 8
        Right (f1,f2) ->
            let f1' = fromSomeTags f1 :: File '[] 'Fastq
                f2' = fromSomeTags f2 :: File '[] 'Fastq
            in bwaAlign output idx (Right (f1', f2')) $
                defaultBWAOpts & bwaCores .~ 8
        )

filterBamSort :: SCATACSeqConfig config
              => SCATACSeq S (Either (File '[] 'Bam) (File '[PairedEnd] 'Bam))
              -> ReaderT config IO ( SCATACSeq S ( Either
                  (File '[NameSorted] 'Bam)
                  (File '[NameSorted, PairedEnd] 'Bam) ))
filterBamSort input = do
    dir <- asks ((<> "/Bam") . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_filt.bam" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . either
        (fmap Left . fun output) (fmap Right . filterBam "./" output)
  where
    fun output x = withTemp (Just "./") $ \f ->
        filterBam "./" f x >>= sortBamByName "./" output

qualityControl :: SCATACSeqConfig config
               => SCATACSeq S ( Either
                   (File '[NameSorted] 'Bam)
                   (File '[NameSorted, PairedEnd] 'Bam) )
               -> ReaderT config IO
                   (SCATACSeq S ( File '[NameSorted, Gzip] 'Bed  -- ^ filtered reads
                             , File '[NameSorted] 'Bam        -- ^ mito reads
                             , File '[] 'Tsv ))               -- ^ Stat
qualityControl input = do
    dir <- asks _scatacseq_output_dir >>= getPath . (<> (asDir "/Bed"))
    anno <- fromJust <$> asks _scatacseq_annotation 
    let output1 = printf "%s/%s_rep%d_filt.bed.gz" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
        output2 = printf "%s/%s_rep%d_stat.txt" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
        mitoReads = printf "%s/%s_rep%d_mito.bam" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
        f fl = do
            header <- getBamHeader fl
            promoters <- (fmap.fmap) (const ()) <$> readPromoters anno
            _ <- runResourceT $ runConduit $ streamBam fl .|
                groupBy ((==) `on` (extractBarcode . queryName)) .|
                mapC (filterReads header promoters) .|
                zipSinks ( zipSinks (tagsOutputSink output1)
                    (concatMapC (snd . fst) .| sinkBam mitoReads header) )
                    (qcFileSink output2)
            return ( location .~ output1 $ emptyFile
                   , location .~ mitoReads $ emptyFile
                   , location .~ output2 $ emptyFile )
    input & replicates.traverse.files %%~ liftIO . f . either
        (^.location) (^.location)
  where
    tagsOutputSink out = concatMapC (fmap fst . filterCell) .|
        concatC .| sinkFileBedGzip out
    qcFileSink out = mapC (showStat . snd) .| unlinesAsciiC .| sinkFile out

filterCell :: (a, Stat) -> Maybe a 
filterCell (x, Stat{..})
    | _uniq_reads >= 1000 && _frip <= 0.8 && _frip >= 0.2 = Just x
    | otherwise = Nothing
{-# INLINE filterCell #-}

-- | Remove duplicated and chrM reads. 
filterReads :: BAMHeader
            -> BEDTree a   -- ^ Promoter
            -> [BAM]
            -> (([BED], [BAM]), Stat)
filterReads hdr promoters bam = runState filterFn $ Stat bc 0 0 0 0
  where
    filterFn = do
        (tags, mito) <- {-rmAbnoramlFragment bam >>= -} rmDup bam >>= rmChrM hdr
        let beds = mapMaybe (fmap changeName . bamToBed hdr) tags
        frip promoters beds
        totalReads beds
        return (beds, mito)
    bc = extractBarcode $ queryName $ head bam
    changeName x = name .~ Just bc $ x
{-# INLINE filterReads #-}
