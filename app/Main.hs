{-# LANGUAGE DeriveGeneric     #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}

module Main where

import           Bio.Pipeline.CallPeaks
import           Bio.Pipeline.Utils
import           Control.Lens                  ((&), (.~))
import           Data.Aeson                    (FromJSON, ToJSON)
import           Data.Default
import           GHC.Generics                  (Generic)
import           Scientific.Workflow
import           Scientific.Workflow.Main      (MainOpts (..), defaultMainOpts,
                                                mainWith)
import           Taiji.Pipeline.SC.ATACSeq        (builder)
import Taiji.Pipeline.SC.ATACSeq.Types

data SCATACSeqOpts = SCATACSeqOpts
    { outputDir :: Directory
    , bwaIndex  :: Maybe FilePath
    , genome    :: Maybe FilePath
    , input     :: FilePath
    , picard    :: Maybe FilePath
    , motifFile :: Maybe FilePath
    , genomeIndex :: Maybe FilePath
    } deriving (Generic)

instance SCATACSeqConfig SCATACSeqOpts where
    _scatacseq_output_dir = outputDir
    _scatacseq_bwa_index = bwaIndex
    _scatacseq_genome_fasta = genome
    _scatacseq_input = input
    _scatacseq_picard = picard
    _scatacseq_callpeak_opts _ = def & mode .~ NoModel (-100) 200
    _scatacseq_genome_index = genomeIndex
    _scatacseq_motif_file = motifFile

instance Default SCATACSeqOpts where
    def = SCATACSeqOpts
        { outputDir = asDir "output"
        , bwaIndex = Nothing
        , genome = Nothing
        , input = "input.yml"
        , picard = Nothing
        , motifFile = Nothing
        , genomeIndex = Nothing
        }

instance FromJSON SCATACSeqOpts
instance ToJSON SCATACSeqOpts

mainWith defaultMainOpts
    { programHeader = "Taiji-ATAC-Seq"
    , workflowConfigType = Just ''SCATACSeqOpts
    } builder
