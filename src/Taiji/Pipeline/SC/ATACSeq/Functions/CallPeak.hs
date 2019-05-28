{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE DataKinds #-}

module Taiji.Pipeline.SC.ATACSeq.Functions.CallPeak
    ( callPeakCluster
    , mkCellClusterBed
    , subSampleClusterBed
    ) where

import           Bio.Pipeline
import Control.Monad (forM)
import Bio.Data.Experiment
import qualified Data.HashSet as S
import Control.Monad.Reader (asks, liftIO)
import Data.Maybe
import Control.Lens
import Bio.Data.Bed
import Conduit
import Data.Conduit.Internal (zipSinks)
import qualified Data.Text as T
import qualified Data.ByteString.Char8 as B
import           Scientific.Workflow
import Shelly hiding (FilePath)

import Taiji.Types
import Taiji.Pipeline.SC.ATACSeq.Types

callPeakCluster :: SCATACSeqConfig config
                => (SCATACSeq S (B.ByteString, File '[Gzip] 'Bed))
                -> WorkflowConfig config
                    (SCATACSeq S (B.ByteString, File '[] 'NarrowPeak))
callPeakCluster input = do
    let idRep = asDir $ "/Peaks/" <> T.unpack (input^.eid) <>
            "_rep" <> show (input^.replicates._1)
    dir <- asks _scatacseq_output_dir >>= getPath . (<> idRep)
    opts <- asks _scatacseq_callpeak_opts
    input & replicates.traverse.files %%~ liftIO . ( \(nm, fl) -> do
        let output = dir ++ "/" ++ B.unpack nm ++ ".narrowPeak" 
        r <- callPeaks output fl Nothing opts
        return (nm, r) )

-- | Extract BEDs for each cluster.
mkCellClusterBed :: SCATACSeqConfig config
                 => SCATACSeq S ( File '[NameSorted, Gzip] 'Bed
                                , [CellCluster] )  -- ^ clusters
                 -> WorkflowConfig config
                    (SCATACSeq S [(B.ByteString, File '[Gzip] 'Bed, Int)])
mkCellClusterBed input = do
    let idRep = asDir $ "/Bed/" <> T.unpack (input^.eid) <>
            "_rep" <> show (input^.replicates._1)
    dir <- asks _scatacseq_output_dir >>= getPath . (<> idRep)
    input & replicates.traverse.files %%~ liftIO . ( \(bed, cs) -> do
        let sinks = sequenceConduits $ flip map cs $ \CellCluster{..} -> do
                let output = dir ++ "/" ++ B.unpack _cluster_name ++ ".bed.gz"
                    cells = S.fromList $ map _cell_id _cluster_member
                    fl = location .~ output $ emptyFile
                (_, depth) <- filterC (f cells) .|
                    zipSinks (sinkFileBedGzip output) lengthC
                return (_cluster_name, fl, depth)
        runResourceT $ runConduit $
            streamBedGzip (bed^.location) .| sinks )
  where
    f :: S.HashSet B.ByteString -> BED -> Bool
    f cells x = fromJust (x^.name) `S.member` cells

-- | Subsampling bed files.
subSampleClusterBed :: SCATACSeqConfig config
                 => (SCATACSeq S [(B.ByteString, File '[Gzip] 'Bed, Int)])
                 -> WorkflowConfig config
                    (SCATACSeq S [(B.ByteString, File '[Gzip] 'Bed)])
subSampleClusterBed input = do
    let idRep = asDir $ "/Bed/" <> T.unpack (input^.eid) <>
            "_rep" <> show (input^.replicates._1) <> "Subsample"
    dir <- asks _scatacseq_output_dir >>= getPath . (<> idRep)
    input & replicates.traverse.files %%~ liftIO . ( \fls -> 
        fmap catMaybes $ forM fls $ \(name, bed, n) -> if n >= depth
            then do
                let output = dir ++ "/" ++ B.unpack name ++ ".bed.gz"
                subSample depth (bed^.location) output
                return $ Just (name, location .~ output $ emptyFile)
            else return Nothing )
  where
    -- decide the sequencing depth
    depth = minimum $ filter (>=lowestDepth) $
        input^..replicates._2.files.folded._3
    lowestDepth = 1000000

subSample :: Int -> FilePath -> FilePath -> IO ()
subSample n input output = shelly $ escaping False $ silently $ 
    bashPipeFail bash_ "zcat"
        [ T.pack input, "|", "shuf", "-n", T.pack $ show n
        , "|", "gzip", "-c", ">", T.pack output ]