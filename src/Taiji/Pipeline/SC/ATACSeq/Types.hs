{-# LANGUAGE DeriveGeneric          #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE UndecidableInstances #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE StandaloneDeriving #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE LambdaCase #-}

module Taiji.Pipeline.SC.ATACSeq.Types
    ( SCATACSeq(..)
    , SCATACSeqConfig(..)
    , qcDir
    , figDir
    , getGenomeFasta
    , getAnnotation
    , getMotif
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
    _scatacseq_assembly :: config -> Maybe String
    _scatacseq_bwa_index :: config -> Maybe FilePath
    _scatacseq_genome_fasta :: config -> Maybe FilePath
    _scatacseq_genome_index :: config -> Maybe FilePath
    _scatacseq_motif_file :: config -> Maybe FilePath
    _scatacseq_callpeak_opts :: config -> CallPeakOpts
    _scatacseq_annotation :: config -> Maybe FilePath
    _scatacseq_tmp_dir :: config -> Maybe FilePath
    _scatacseq_cluster_resolution :: config -> Maybe Double
    _scatacseq_blacklist :: config -> Maybe FilePath
    _scatacseq_te_cutoff :: config -> Double

qcDir :: SCATACSeqConfig config => ReaderT config IO FilePath
qcDir = asks _scatacseq_output_dir >>= getPath . (<> "/QC/")

figDir :: SCATACSeqConfig config => ReaderT config IO FilePath
figDir = asks _scatacseq_output_dir >>= getPath . (<> "/Figure/")

getGenomeFasta :: SCATACSeqConfig config => ReaderT config IO FilePath
getGenomeFasta = asks _scatacseq_genome_fasta >>= \case
    Nothing -> error "genome fasta is missing"
    Just fasta -> do
       exist <- liftIO $ shelly $ test_f $ fromText $ T.pack fasta 
       if exist
           then return fasta
           else asks _scatacseq_assembly >>= \case
               Nothing -> error "genome fasta is missing"
               Just assembly -> do
                   liftIO $ fetchGenome fasta assembly
                   return fasta

getAnnotation :: SCATACSeqConfig config => ReaderT config IO FilePath
getAnnotation = asks _scatacseq_annotation >>= \case
    Nothing -> error "annotation is missing"
    Just anno -> do
       exist <- liftIO $ shelly $ test_f $ fromText $ T.pack anno
       if exist
           then return anno
           else asks _scatacseq_assembly >>= \case
               Nothing -> error "annotation is missing"
               Just assembly -> do
                   liftIO $ fetchAnnotation anno assembly
                   return anno

getMotif :: SCATACSeqConfig config => ReaderT config IO FilePath
getMotif = asks _scatacseq_motif_file >>= \case
    Nothing -> error "motif file is missing"
    Just motif -> do
       exist <- liftIO $ shelly $ test_f $ fromText $ T.pack motif
       if exist
           then return motif
           else asks _scatacseq_assembly >>= \case
               Nothing -> error "motif file is missing"
               Just assembly -> do
                   liftIO $ fetchMotif motif assembly
                   return motif

getGenomeIndex :: SCATACSeqConfig config => ReaderT config IO FilePath
getGenomeIndex = do
    seqIndex <- asks ( fromMaybe (error "Genome index file was not specified!") .
        _scatacseq_genome_index )
    fileExist <- shelly $ test_f $ fromText $ T.pack seqIndex
    unless fileExist $ do
        genome <- getGenomeFasta
        shelly $ mkdir_p $ fromText $ T.pack $ takeDirectory seqIndex
        liftIO $ mkIndex [genome] seqIndex
    return seqIndex
{-# INLINE getGenomeIndex #-}