{-# LANGUAGE DeriveGeneric          #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE UndecidableInstances #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE StandaloneDeriving #-}
{-# LANGUAGE OverloadedStrings #-}

module Taiji.Pipeline.SC.ATACSeq.Types
    ( SCATACSeq(..)
    , SCATACSeqConfig(..)
    , qcDir
    , figDir
    , getGenomeIndex
    ) where

import           Data.Binary (Binary(..))
import Bio.Data.Experiment.Types
import Bio.Data.Experiment.Replicate
import GHC.Generics (Generic)
import qualified Data.Text as T
import           Bio.Seq.IO
import           System.FilePath               (takeDirectory)
import           Shelly                        (fromText, mkdir_p, shelly,
                                                test_f)

import Taiji.Prelude

newtype SCATACSeq container file = SCATACSeq (CommonFields container file)
     deriving (Generic, Experiment)

deriving instance Show (container (Replicate file)) => Show (SCATACSeq container file)

instance Binary (container (Replicate file)) =>
    Binary (SCATACSeq container file)

class SCATACSeqConfig config where
    _scatacseq_output_dir :: config -> Directory
    _scatacseq_input :: config -> FilePath
    _scatacseq_bwa_index :: config -> Maybe FilePath
    _scatacseq_genome_fasta :: config -> Maybe FilePath
    _scatacseq_genome_index :: config -> Maybe FilePath
    _scatacseq_motif_file :: config -> Maybe FilePath
    _scatacseq_callpeak_opts :: config -> CallPeakOpts
    _scatacseq_annotation :: config -> Maybe FilePath
    _scatacseq_temp_dir :: config -> Maybe FilePath
    _scatacseq_cluster_resolution :: config -> Maybe Double
    _scatacseq_marker_gene_list :: config -> Maybe FilePath
    _scatacseq_blacklist :: config -> Maybe FilePath

qcDir :: SCATACSeqConfig config => ReaderT config IO FilePath
qcDir = asks _scatacseq_output_dir >>= getPath . (<> "/QC/")

figDir :: SCATACSeqConfig config => ReaderT config IO FilePath
figDir = asks _scatacseq_output_dir >>= getPath . (<> "/Figure/")

getGenomeIndex :: SCATACSeqConfig config => ReaderT config IO FilePath
getGenomeIndex = do
    seqIndex <- asks ( fromMaybe (error "Genome index file was not specified!") .
        _scatacseq_genome_index )
    genome <- asks ( fromMaybe (error "Genome fasta file was not specified!") .
        _scatacseq_genome_fasta )
    shelly $ do
        fileExist <- test_f $ fromText $ T.pack seqIndex
        unless fileExist $ do
            mkdir_p $ fromText $ T.pack $ takeDirectory seqIndex
            liftIO $ mkIndex [genome] seqIndex
    return seqIndex
{-# INLINE getGenomeIndex #-}

