{-# LANGUAGE DeriveGeneric          #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE UndecidableInstances #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
module Taiji.Pipeline.SC.ATACSeq.Types
    ( SCATACSeq(..)
    , SCATACSeqConfig(..)
    ) where

import           Data.Serialize                (Serialize (..))
import           Data.Aeson
import Bio.Data.Experiment.Types
import Bio.Data.Experiment.Replicate
import Bio.Pipeline.Utils (Directory)
import Bio.Pipeline.CallPeaks (CallPeakOpts)
import GHC.Generics (Generic)

newtype SCATACSeq container file = SCATACSeq (CommonFields container file)
     deriving (Generic, Experiment)

instance FromJSON (container (Replicate file)) =>
    FromJSON (SCATACSeq container file)

instance ToJSON (container (Replicate file)) =>
    ToJSON (SCATACSeq container file)

instance Serialize (container (Replicate file)) =>
    Serialize (SCATACSeq container file)

class SCATACSeqConfig config where
    _scatacseq_output_dir :: config -> Directory
    _scatacseq_input :: config -> FilePath
    _scatacseq_bwa_index :: config -> Maybe FilePath
    _scatacseq_genome_fasta :: config -> Maybe FilePath
    _scatacseq_genome_index :: config -> Maybe FilePath
    _scatacseq_motif_file :: config -> Maybe FilePath
    _scatacseq_callpeak_opts :: config -> CallPeakOpts
