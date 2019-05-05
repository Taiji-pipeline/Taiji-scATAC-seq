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
    { output_dir :: Directory
    , bwa_index  :: Maybe FilePath
    , genome :: Maybe FilePath
    , input     :: FilePath
    , motif_file :: Maybe FilePath
    , genome_index :: Maybe FilePath
    , annotation :: Maybe FilePath
    } deriving (Generic)

instance SCATACSeqConfig SCATACSeqOpts where
    _scatacseq_output_dir = output_dir
    _scatacseq_bwa_index = bwa_index
    _scatacseq_genome_fasta = genome
    _scatacseq_input = input
    _scatacseq_callpeak_opts _ = def & mode .~ NoModel (-100) 200
    _scatacseq_genome_index = genome_index
    _scatacseq_motif_file = motif_file
    _scatacseq_annotation = annotation

instance Default SCATACSeqOpts where
    def = SCATACSeqOpts
        { output_dir = asDir "output"
        , bwa_index = Nothing
        , genome_index = Nothing
        , input = "input.yml"
        , motif_file = Nothing
        , genome = Nothing
        , annotation = Nothing
        }

instance FromJSON SCATACSeqOpts
instance ToJSON SCATACSeqOpts

mainWith defaultMainOpts
    { programHeader = "Taiji-ATAC-Seq"
    , workflowConfigType = Just ''SCATACSeqOpts
    } builder
