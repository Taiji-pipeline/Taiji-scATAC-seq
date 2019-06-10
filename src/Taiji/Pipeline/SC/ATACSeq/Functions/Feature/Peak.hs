{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE DataKinds #-}

module Taiji.Pipeline.SC.ATACSeq.Functions.Feature.Peak
    ( findPeaks
    , mergePeaks
    , mkCellClusterBed
    , subSampleClusterBed
    ) where

import           Bio.Pipeline
import qualified Data.HashSet as S
import Bio.Data.Bed
import Data.Conduit.Internal (zipSinks)
import qualified Data.Text as T
import qualified Data.ByteString.Char8 as B
import Shelly hiding (FilePath)

import Taiji.Prelude
import Taiji.Pipeline.SC.ATACSeq.Types

{-
callPeakBulk :: SCATACSeqConfig config
             => SCATACSeq S (File tags 'Bed)
             -> ReaderT config IO (SCATACSeq S (File '[] 'NarrowPeak))
callPeakBulk input = do
    dir <- asks _scatacseq_output_dir >>= getPath . (<> (asDir "/Peaks"))
    let output = printf "%s/%s_rep%d_filt.bed.gz" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    opts <- asks _scatacseq_callpeak_opts
    input & replicates.traverse.files %%~ liftIO . ( \(nm, fl) -> do
        let output = dir ++ "/" ++ B.unpack nm ++ ".narrowPeak" 
        r <- callPeaks output fl Nothing opts
        return (nm, r) )
        -}

-- | Call Peaks
findPeaks :: SCATACSeqConfig config
                => (B.ByteString, File '[Gzip] 'Bed)
                -> ReaderT config IO (B.ByteString, File '[] 'NarrowPeak)
findPeaks (cName, bedFl) = do
    dir <- asks _scatacseq_output_dir >>= getPath . (<> "/Peaks/Cluster/")
    opts <- asks _scatacseq_callpeak_opts
    let output = dir ++ B.unpack cName ++ ".narrowPeak" 
    r <- liftIO $ callPeaks output bedFl Nothing opts
    return (cName, r)

mergePeaks :: SCATACSeqConfig config
           => [(B.ByteString, File '[] 'NarrowPeak)]
           -> ReaderT config IO (Maybe (File '[Gzip] 'NarrowPeak))
mergePeaks [] = return Nothing
mergePeaks input = do
    dir <- asks _scatacseq_output_dir >>= getPath . (<> "/Peaks/")
    let output = dir <> "merged.narrowPeak.gz" 
    liftIO $ withTemp Nothing $ \tmp -> do
        shelly $ escaping False $ bashPipeFail bash_ "cat" $
            map (T.pack . (^._2.location)) input ++
            [ "|", "sort", "-k1,1", "-k2,2n", "-k3,3n", ">", T.pack tmp ]
        runResourceT $ runConduit $ streamBed tmp .| 
            mergeSortedBedWith getBestPeak .| mapC resize .| sinkFileBedGzip output
    return $ Just $ location .~ output $ emptyFile
  where
    getBestPeak = maximumBy $ comparing (fromJust . (^.npPvalue))
    resize pk = chromStart .~ max 0 (summit - 250) $
        chromEnd .~ summit + 250 $ pk
      where
        summit = pk^.chromStart + fromJust (pk^.npPeak)
    
-- | Extract BEDs for each cluster.
mkCellClusterBed :: SCATACSeqConfig config
                 => SCATACSeq S ( File '[NameSorted, Gzip] 'Bed
                                , [CellCluster] )  -- ^ clusters
                 -> ReaderT config IO
                    (SCATACSeq S [(B.ByteString, File '[Gzip] 'Bed, Int)])
mkCellClusterBed input = do
    let idRep = asDir $ "/Bed/" <> T.unpack (input^.eid) <>
            "_rep" <> show (input^.replicates._1)
    dir <- asks _scatacseq_output_dir >>= getPath . (<> idRep)
    input & replicates.traverse.files %%~ liftIO . ( \(bed, cs) -> do
        let sinks = sequenceConduits $ flip map cs $ \CellCluster{..} -> do
                let output = dir ++ "/" ++ B.unpack _cluster_name ++ ".bed.gz"
                    cells = S.fromList $ map _cell_barcode _cluster_member
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
                 -> ReaderT config IO
                    (SCATACSeq S [(B.ByteString, File '[Gzip] 'Bed)])
subSampleClusterBed input = do
    let idRep = asDir $ "/Bed/" <> T.unpack (input^.eid) <>
            "_rep" <> show (input^.replicates._1) <> "_Subsample"
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