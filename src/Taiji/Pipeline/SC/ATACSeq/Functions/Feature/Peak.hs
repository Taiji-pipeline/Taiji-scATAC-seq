{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE GADTs #-}

module Taiji.Pipeline.SC.ATACSeq.Functions.Feature.Peak
    ( mkPeakMat
    , findPeaks
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
import Data.Singletons.Prelude (Elem)
import Shelly hiding (FilePath)

import Taiji.Prelude
import Taiji.Pipeline.SC.ATACSeq.Types
import Taiji.Pipeline.SC.ATACSeq.Functions.Utils

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

-- | Make the read count matrix.
mkPeakMat :: (Elem 'Gzip tags ~ 'True, SCATACSeqConfig config)
          => SCATACSeq S (File tags 'Bed, File '[Gzip] 'NarrowPeak, Int)
          -> ReaderT config IO (SCATACSeq S (File tags 'Other))
mkPeakMat input = do
    dir <- asks ((<> "/Feature/") . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_peak.txt.gz" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . (\(tagFl, regionFl, nCell) -> do
        regions <- runResourceT $ runConduit $
            streamBedGzip (regionFl^.location) .| sinkList :: IO [BED3]
        runResourceT $ runConduit $ streamBedGzip (tagFl^.location) .|
            groupCells .| mkFeatMat nCell regions .| sinkFile output
        return $ emptyFile & location .~ output )

-- | Call Peaks for aggregated files
findPeaks :: SCATACSeqConfig config
          => FilePath
          -> (B.ByteString, File '[Gzip] 'Bed)
          -> ReaderT config IO (B.ByteString, File '[] 'NarrowPeak)
findPeaks prefix (cName, bedFl) = do
    dir <- asks _scatacseq_output_dir >>= getPath . (<> asDir prefix)
    opts <- asks _scatacseq_callpeak_opts
    let output = dir ++ "/" ++ B.unpack cName ++ ".narrowPeak" 
    r <- liftIO $ callPeaks output bedFl Nothing opts
    return (cName, r)

mergePeaks :: SCATACSeqConfig config
           => FilePath
           -> [(B.ByteString, File '[] 'NarrowPeak)]
           -> ReaderT config IO (Maybe (File '[Gzip] 'NarrowPeak))
mergePeaks _ [] = return Nothing
mergePeaks prefix input = do
    dir <- asks _scatacseq_output_dir >>= getPath . (<> asDir prefix)
    let output = dir <> "/merged.narrowPeak.gz" 
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