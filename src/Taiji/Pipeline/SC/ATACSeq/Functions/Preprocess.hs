{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE DataKinds #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Preprocess
    ( readInput
    , downloadData
    , getFastq
    , getSortedBed
    ) where

import Bio.Data.Experiment
import Control.Lens
import           Data.Bifunctor                (bimap)
import Bio.Pipeline
import Data.Either (rights, lefts)
import Bio.Data.Experiment.Parser
import Scientific.Workflow
import Control.Monad.Reader (asks, liftIO)

import Taiji.Pipeline.SC.ATACSeq.Types

type RAWInput = SCATACSeq N [Either SomeFile (SomeFile, SomeFile)]

readInput :: SCATACSeqConfig config
          => () -> WorkflowConfig config [RAWInput]
readInput _ = do
    input <- asks _scatacseq_input
    liftIO $ simpleInputReader input "scATAC-seq" SCATACSeq

downloadData :: SCATACSeqConfig config
             => [RAWInput]
             -> WorkflowConfig config [RAWInput]
downloadData input = do
    dir <- asks _scatacseq_output_dir >>= getPath . (<> (asDir "/Download"))
    liftIO $ input & traverse.replicates.traverse.files.traverse %%~
        downloadFiles dir

getFastq :: [RAWInput]
         -> [ SCATACSeq S (SomeTags 'Fastq, SomeTags 'Fastq) ]
getFastq input = concatMap split $ concatMap split $
    input & mapped.replicates.mapped.files %~ f
  where
    f fls = map (bimap castFile castFile) $
        filter (\(x,y) -> getFileType x == Fastq && getFileType y == Fastq) $
        rights fls

getSortedBed :: [RAWInput]
             -> [ SCATACSeq S (File '[NameSorted, Gzip] 'Bed) ]
getSortedBed input = concatMap split $ concatMap split $
    input & mapped.replicates.mapped.files %~ f
  where
    f fls = map fromSomeFile $
        filter (\x -> getFileType x == Bed && x `hasTag` NameSorted &&
            x `hasTag` Gzip) $ lefts fls