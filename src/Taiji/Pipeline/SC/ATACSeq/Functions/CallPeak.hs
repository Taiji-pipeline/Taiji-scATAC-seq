{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE DataKinds #-}

module Taiji.Pipeline.SC.ATACSeq.Functions.CallPeak
    ( callPeakCluster
    , mkCellClusterBed
    ) where

import           Bio.Pipeline
import Bio.Data.Experiment
import qualified Data.HashSet as S
import Control.Monad.Reader (asks, liftIO)
import Data.Maybe
import Control.Lens
import Bio.Data.Bed
import Conduit
import qualified Data.Text as T
import qualified Data.ByteString.Char8 as B
import           Scientific.Workflow

import Taiji.Pipeline.SC.ATACSeq.Functions.Utils
import Taiji.Pipeline.SC.ATACSeq.Types

callPeakCluster :: SCATACSeqConfig config
                => (SCATACSeq S [(B.ByteString, File '[NameSorted, Gzip] 'Bed)])
                -> WorkflowConfig config
                    (SCATACSeq S [(B.ByteString, File '[] 'NarrowPeak)])
callPeakCluster input = do
    let idRep = asDir $ "/Peaks/" <> T.unpack (input^.eid) <>
            "_rep" <> show (input^.replicates._1)
    dir <- asks _scatacseq_output_dir >>= getPath . (<> idRep)
    opts <- asks _scatacseq_callpeak_opts
    input & replicates.traverse.files.traverse %%~ liftIO . ( \(nm, fl) -> do
        let output = dir ++ "/" ++ B.unpack nm ++ ".narrowPeak" 
        r <- callPeaks output fl Nothing opts
        return (nm, r) )

-- | Extract BEDs for each cluster.
mkCellClusterBed :: SCATACSeqConfig config
                 => SCATACSeq S ( File '[NameSorted, Gzip] 'Bed
                                , File '[] 'Other )  -- ^ clusters
                 -> WorkflowConfig config
                    (SCATACSeq S [(B.ByteString, File '[NameSorted, Gzip] 'Bed)])
mkCellClusterBed input = do
    let idRep = asDir $ "/Bed/" <> T.unpack (input^.eid) <>
            "_rep" <> show (input^.replicates._1)
    dir <- asks _scatacseq_output_dir >>= getPath . (<> idRep)
    input & replicates.traverse.files %%~ liftIO . ( \(bed, cluster) -> do
        cs <- readCellCluster $ cluster^.location
        let sinks = sequenceConduits $ flip map cs $ \CellCluster{..} ->
                let output = dir ++ "/" ++ B.unpack _cluster_name ++ ".bed.gz"
                    cells = S.fromList _cluster_member
                    fl = location .~ output $ emptyFile
                in filterC (f cells) .|
                      sinkFileBedGzip output >> return (_cluster_name, fl)
        runResourceT $ runConduit $
            streamBedGzip (bed^.location) .| sinks )
  where
    f :: S.HashSet B.ByteString -> BED -> Bool
    f cells x = fromJust (x^.name) `S.member` cells
        
