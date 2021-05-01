{-# LANGUAGE DeriveGeneric          #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE UndecidableInstances #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE StandaloneDeriving #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE DeriveAnyClass #-}

module Taiji.Pipeline.SC.ATACSeq.Types
    ( SCATACSeq(..)
    , SCATACSeqConfig(..)
    , Stat(..)
    , TagsAligned
    , qcDir
    , figDir
    , getGenomeFasta
    , getAnnotation
    , getMotif
    , mkGenomeIndex
    , getCallPeakOpt
    , getQCFunction
    ) where

import Numeric.Natural (Natural)
import           Data.Binary (Binary(..))
import Bio.Data.Experiment.Types
import qualified Data.ByteString.Char8 as B
import Bio.Data.Experiment.Replicate
import GHC.Generics (Generic)
import qualified Data.Map.Strict as M
import qualified Data.Text as T
import Control.Exception (catch, SomeException(..))
import           Bio.Seq.IO
import           System.FilePath               (takeDirectory)
import           Shelly                        (fromText, mkdir_p, shelly,
                                                test_f)
import Control.DeepSeq (NFData)

import Taiji.Prelude

newtype SCATACSeq container file = SCATACSeq (CommonFields container file)
    deriving Generic
    deriving newtype Experiment

deriving instance Show (container (Replicate file)) => Show (SCATACSeq container file)

instance Binary (container (Replicate file)) =>
    Binary (SCATACSeq container file)

type TagsAligned =
    ( Either (File '[NameSorted, Gzip] 'Bed) (File '[NameSorted, PairedEnd, Gzip] 'Bed)
    , File '[] 'Tsv )

class SCATACSeqConfig config where
    _scatacseq_output_dir :: config -> Directory
    _scatacseq_input :: config -> FilePath
    _scatacseq_batch_info :: config -> Maybe FilePath
    _scatacseq_assembly :: config -> Maybe String
    _scatacseq_bwa_index :: config -> Maybe FilePath
    _scatacseq_genome_fasta :: config -> Maybe FilePath
    _scatacseq_genome_index :: config -> Maybe FilePath
    _scatacseq_motif_file :: config -> Maybe FilePath
    _scatacseq_callpeak_opts :: config -> CallPeakOpts
    _scatacseq_annotation :: config -> Maybe FilePath
    _scatacseq_tmp_dir :: config -> Maybe FilePath
    _scatacseq_cell_barcode_length :: config -> Maybe Natural
    _scatacseq_cluster_resolution_list :: config -> [Double]
    _scatacseq_cluster_resolution :: config -> Maybe Double
    _scatacseq_cluster_exclude :: config -> [T.Text]
    _scatacseq_do_subclustering :: config -> Bool
    _scatacseq_subcluster_resolution :: config -> Maybe (M.Map T.Text Double)
    _scatacseq_cluster_optimizer :: config -> Optimizer
    _scatacseq_blacklist :: config -> Maybe FilePath
    _scatacseq_te_cutoff :: config -> Double
    _scatacseq_minimal_fragment :: config -> Natural
    _scatacseq_doublet_score_cutoff :: config -> Double
    _scatacseq_cluster_by_window :: config -> Bool
    _scatacseq_window_size :: config -> Natural

data Stat = Stat
    { _barcode :: B.ByteString
    , _dup_rate :: Maybe Double
    , _mito_rate :: Maybe Double
    , _te :: Double
    , _uniq_reads :: Int
    , _doublet_score :: Maybe Double } deriving (Generic, NFData)

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

mkGenomeIndex :: SCATACSeqConfig config => ReaderT config IO ()
mkGenomeIndex = do
    seqIndex <- asks ( fromMaybe (error "Genome index file was not specified!") .
        _scatacseq_genome_index )
    exist <- liftIO $ catch (withGenome seqIndex $ const $ return True) $
        \(SomeException _) -> return False
    unless exist $ do
        liftIO $ putStrLn "Generating genome index ..."
        genome <- getGenomeFasta
        shelly $ mkdir_p $ fromText $ T.pack $ takeDirectory seqIndex
        liftIO $ mkIndex [genome] seqIndex
{-# INLINE mkGenomeIndex #-}

getCallPeakOpt :: SCATACSeqConfig config => ReaderT config IO CallPeakOpts
getCallPeakOpt = do
    opt <- asks _scatacseq_callpeak_opts 
    s <- case opt^.gSize of
        Nothing -> do
            idx <- asks ( fromMaybe (error "Genome index file was not specified!") .
                _scatacseq_genome_index )
            s <- liftIO $ fmap fromIntegral $ withGenome idx $
                return . foldl1' (+) . map snd . getChrSizes
            return $ Just $ show (truncate $ 0.9 * (s :: Double) :: Int)
        x -> return x
    return $ gSize .~ s $ opt
{-# INLINE getCallPeakOpt #-}

getQCFunction :: SCATACSeqConfig config => ReaderT config IO (Stat -> Bool)
getQCFunction = do
    teCutoff <- asks _scatacseq_te_cutoff
    fragmentCutoff <- fromIntegral <$> asks _scatacseq_minimal_fragment
    doubletCutoff <- asks _scatacseq_doublet_score_cutoff
    return $ \x -> _te x >= teCutoff && _uniq_reads x >= fragmentCutoff && fromMaybe 0 (_doublet_score x) <= doubletCutoff
{-# INLINE getQCFunction #-}