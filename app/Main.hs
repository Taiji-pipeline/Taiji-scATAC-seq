{-# LANGUAGE DeriveGeneric     #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}

module Main where

import           Bio.Pipeline.CallPeaks
import           Bio.Pipeline.Utils
import           Control.Lens                  ((&), (.~))
import           Data.Aeson                    (FromJSON)
import Data.Binary (Binary)
import           GHC.Generics                  (Generic)
import Data.Default

import           Control.Workflow
import qualified Control.Workflow.Coordinator.Drmaa as D
import Control.Workflow.Main
import Data.Proxy (Proxy(..))

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

instance Binary SCATACSeqOpts
instance FromJSON SCATACSeqOpts

instance SCATACSeqConfig SCATACSeqOpts where
    _scatacseq_output_dir = output_dir
    _scatacseq_bwa_index = fmap (++ "/genome.fa") . bwa_index
    _scatacseq_genome_fasta = genome
    _scatacseq_input = input
    _scatacseq_callpeak_opts _ = def & mode .~ NoModel (-100) 200
                                   & cutoff .~ QValue 0.05
                                   & callSummits .~ True
    _scatacseq_genome_index = genome_index
    _scatacseq_motif_file = motif_file
    _scatacseq_annotation = annotation

decodeDrmaa :: String -> Int -> FilePath -> IO D.DrmaaConfig
decodeDrmaa ip port _ = D.getDefaultDrmaaConfig
    ["remote", "--ip", ip, "--port", show port]

build "wf" [t| SciFlow SCATACSeqOpts |] builder

main :: IO ()
main = defaultMain "" cmd wf
  where
    cmd = [ runParser decodeDrmaa
          , viewParser
          , deleteParser
          , remoteParser (Proxy :: Proxy D.Drmaa) ]

