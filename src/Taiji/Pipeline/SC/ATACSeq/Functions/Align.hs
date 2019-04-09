{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Align
    ( tagAlign
    , filterBamSort
    , rmDuplicates
    ) where

import Bio.Data.Experiment
import           Bio.HTS
import qualified Data.ByteString.Char8 as B
import Conduit
import Control.Lens
import Bio.Pipeline
import Data.Either (fromRight)
import qualified Data.Text as T
import qualified Data.Map.Strict as M
import Data.Maybe
import Text.Printf (printf)
import Bio.Data.Experiment.Parser
import Scientific.Workflow
import Control.Monad.Reader (asks, liftIO)

import Taiji.Pipeline.SC.ATACSeq.Types

type RAWInput = SCATACSeq N [Either SomeFile (SomeFile, SomeFile)]

tagAlign :: SCATACSeqConfig config
         => SCATACSeq S (SomeTags 'Fastq, SomeTags 'Fastq)
         -> WorkflowConfig config ( SCATACSeq S (File '[PairedEnd] 'Bam) )
tagAlign input = do
    dir <- asks ((<> "/Bam") . _scatacseq_output_dir) >>= getPath
    idx <- asks (fromJust . _scatacseq_bwa_index)
    let output = printf "%s/%s_rep%d.bam" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . ( \(f1,f2) ->
        let f1' = fromSomeTags f1 :: File '[] 'Fastq
            f2' = fromSomeTags f2 :: File '[] 'Fastq
        in fromRight undefined <$> bwaAlign output idx (Right (f1', f2'))
            (defaultBWAOpts & bwaCores .~ 16)
        )

filterBamSort :: SCATACSeqConfig config
              => SCATACSeq S (File '[PairedEnd] 'Bam)
              -> WorkflowConfig config
                  (SCATACSeq S (File '[NameSorted, PairedEnd] 'Bam))
filterBamSort input = do
    dir <- asks ((<> "/Bam") . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_filt.bam" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . filterBam "./" output

rmDuplicates :: SCATACSeqConfig config
             => SCATACSeq S (File '[NameSorted, PairedEnd] 'Bam)
             -> WorkflowConfig config (SCATACSeq S (File '[] 'Bam))
rmDuplicates input = do
    dir <- asks _scatacseq_output_dir >>= getPath . (<> (asDir "/Bam"))
    let output = printf "%s/%s_rep%d_filt_dedup.bam" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
        f fl = do
            header <- getBamHeader fl
            runResourceT $ runConduit $ streamBam fl .| getReadGroup .|
                concatMapC rmDup .| sinkBam output header
            return $ location .~ output $ emptyFile
    input & replicates.traverse.files %%~ liftIO . f . (^.location)

getReadGroup:: Monad m => ConduitT BAM [BAM] m ()
getReadGroup = await >>= \case
    Nothing -> return ()
    Just b -> go (getBarcode b) [b]
  where
    go bc acc = await >>= \case 
        Nothing -> yield acc
        Just b -> do
            let bc' = getBarcode b
            if bc == bc'
                then go bc $ b : acc
                else yield acc >> go bc' [b]
    getBarcode bam = Just $ head $ B.split ':' $ queryName bam

-- | Remove duplicates for reads originated from a single cell.
rmDup :: [BAM] -> [BAM]
rmDup bam = map snd $ M.elems $ M.fromListWith collapse $
    map (\b -> (getKey b, (getSc b, b))) bam
  where
    collapse x y | fst x > fst y = x
                 | otherwise = y
    getKey b | not (hasMultiSegments flg) || isNextUnmapped flg = singleKey
             | otherwise = pairKey
      where
        (singleKey, pairKey) = makeKey (const Nothing) b
        flg = flag b
    getSc b = case queryAuxData ('m', 's') b of
        Just (AuxInt x) -> x + fromJust (sumQual 15 b)
        _ -> fromJust $ sumQual 15 b
{-# INLINE rmDup #-}
