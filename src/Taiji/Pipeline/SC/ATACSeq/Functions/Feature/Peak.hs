{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE FlexibleContexts #-}

module Taiji.Pipeline.SC.ATACSeq.Functions.Feature.Peak
    ( mkPeakMat
    , findPeaks
    , mergePeaks
    , getPeakEnrichment
    , mkCellClusterBed
    , computePeakRAS
    ) where

import           Bio.Pipeline
import qualified Data.HashSet as S
import Bio.Data.Bed
import Data.Conduit.Internal (zipSinks, zipSources)
import qualified Data.Text as T
import qualified Data.ByteString.Char8 as B
import qualified Data.Matrix as Mat
import Data.Singletons.Prelude (Elem, SingI)
import Shelly hiding (FilePath)
import qualified Data.Vector as V

import Statistics.Distribution (quantile)
import Statistics.Distribution.ChiSquared (chiSquared)

import Taiji.Prelude
import Taiji.Pipeline.SC.ATACSeq.Types
import Taiji.Pipeline.SC.ATACSeq.Functions.Utils
import qualified Taiji.Utils.DataFrame as DF

-- | Make the read count matrix.
mkPeakMat :: (Elem 'Gzip tags ~ 'True, SCATACSeqConfig config)
          => FilePath
          -> SCATACSeq S (File tags 'Bed, File '[Gzip] 'NarrowPeak, Int)
          -> ReaderT config IO (SCATACSeq S (File '[Gzip] 'Other))
mkPeakMat prefix input = do
    dir <- asks ((<> asDir prefix) . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_peak.mat.gz" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . (\(tagFl, regionFl, nCell) -> do
        regions <- runResourceT $ runConduit $
            streamBedGzip (regionFl^.location) .| sinkList :: IO [BED3]
        runResourceT $ runConduit $ streamBedGzip (tagFl^.location) .|
            groupCells .| mkFeatMat nCell (map return regions) .| sinkFile output
        return $ emptyFile & location .~ output )

-- | Call Peaks for aggregated files
findPeaks :: (SingI tags, SCATACSeqConfig config)
          => FilePath
          -> CallPeakOpts
          -> (B.ByteString, File tags 'Bed)
          -> ReaderT config IO (B.ByteString, File '[] 'NarrowPeak)
findPeaks prefix opts (cName, bedFl) = do
    dir <- asks _scatacseq_output_dir >>= getPath . (<> asDir prefix)
    let output = dir ++ "/" ++ B.unpack cName ++ ".narrowPeak" 
    asks _scatacseq_blacklist >>= \case
        Nothing -> liftIO $ do
            r <- callPeaks output bedFl Nothing opts
            return (cName, r)
        Just blacklist -> liftIO $ withTemp (Just "./") $ \tmp -> do
            _ <- callPeaks tmp bedFl Nothing opts
            blackRegions <- readBed blacklist :: IO [BED3]
            let bedTree = bedToTree const $ map (\x -> (x, ())) blackRegions
            peaks <- readBed tmp :: IO [NarrowPeak]
            writeBed output $ filter (not . isIntersected bedTree) peaks
            return (cName, location .~ output $ emptyFile)

-- | Merge peaks
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
        let source = streamBed tmp .| mergeSortedBedWith getBestPeak .| mapC resize 
        runResourceT $ runConduit $ zipSources (iterateC succ (0 :: Int)) source .|
            mapC (\(i, p) -> name .~ Just ("p" <> B.pack (show i)) $ p) .|
            sinkFileBedGzip output
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

-- | Get signal enrichment for each peak
getPeakEnrichment :: SCATACSeqConfig config
                  => ( [(B.ByteString, File '[Gzip] 'NarrowPeak)]
                     , Maybe (File '[Gzip] 'NarrowPeak) )
                  -> ReaderT config IO (File '[] 'Tsv, File '[] 'Tsv)
getPeakEnrichment (peaks, Just refPeak) = do
    dir <- asks _scatacseq_output_dir >>= getPath . (<> "/Feature/Peak/")
    let output1 = dir <> "peak_signal.tsv"
        output2 = dir <> "peak_pvalue.tsv"
    liftIO $ do
        list <- runResourceT $ runConduit $
            streamBedGzip (refPeak^.location) .| sinkList
        let col = map toString list
        (row, val) <- fmap unzip $ forM peaks $ \(nm, p) -> do
            peak <- runResourceT $ runConduit $ streamBedGzip (p^.location) .| sinkList
            return (T.pack $ B.unpack nm, getValue peak list)
        let (signal, pval) = DF.unzip $ DF.transpose $ DF.mkDataFrame row col val
        DF.writeTable output1 (T.pack . show) signal
        DF.writeTable output2 (T.pack . show) pval
        return (location .~ output1 $ emptyFile, location .~ output2 $ emptyFile)
  where
    toString x = T.pack (B.unpack $ x^.chrom) <> ":" <>
        T.pack (show $ x^.chromStart) <> "-" <> T.pack (show $ x^.chromEnd)
    getValue :: [NarrowPeak]  -- ^ query
            -> [BED3]  -- ^ Peak list
            -> [(Double, Double)]
    getValue query peakList = runIdentity $ runConduit $
        yieldMany peakList .| intersectBedWith f query .| sinkList
      where
        f _ [] = (0, 0)
        f _ [x] = (x^.npSignal, quantile (chiSquared 1) $ 10**(negate $ fromJust $ x^.npPvalue))
        f _ _ = undefined
getPeakEnrichment _ = undefined

{-
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
        fmap catMaybes $ forM fls $ \(nm, bed, n) -> if n >= depth
            then do
                let output = dir ++ "/" ++ B.unpack nm ++ ".bed.gz"
                subSample depth (bed^.location) output
                return $ Just (nm, location .~ output $ emptyFile)
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

clusterCorrelation :: SCATACSeqConfig config
                   => ( [(B.ByteString, File '[] 'NarrowPeak)]
                      , Maybe (File '[Gzip] 'NarrowPeak) )
                   -> ReaderT config IO ()
clusterCorrelation (input, Just refPeak) = do
    dir <- asks _scatacseq_output_dir >>= getPath . (<> "/Figure/")
    let output = dir <> "cluster_correlation.html"
    liftIO $ do
        ref <- runResourceT $ runConduit $
            streamBedGzip (refPeak^.location) .| sinkList
        peaks <- mapM (readBed . (^.location)) peakFls
        savePlots output [] [ heatmap $ DF.orderDataFrame id $
            DF.mkDataFrame names' names' $ peakCor peaks ref ]
  where
    names' = map (T.pack . B.unpack) names
    (names, peakFls) = unzip input
clusterCorrelation _ = return ()

peakCor :: [[NarrowPeak]] 
        -> [BED3]  -- ^ Peak list
        -> [[Double]]
peakCor peaks refPeak = MU.toLists $ MU.generate (n, n) $ \(i,j) ->
    jaccardIndex (values V.! i) (values V.! j)
  where
    n = V.length values
    values = V.fromList $ map (getValue refPeak) peaks
    getValue peakList query = runST $ runConduit $ yieldMany peakList .|
        intersectBedWith f query .| sinkVector
      where
        f _ [] = False
        f _ _ = True
    jaccardIndex x y = fromIntegral (U.length $ U.filter id $ U.zipWith (&&) x y) /
        fromIntegral (U.length $ U.filter id $ U.zipWith (||) x y)
-}

-- | Compute the relative accessibility score.
computePeakRAS :: SCATACSeqConfig config
               => FilePath
               -> ( Maybe (File '[Gzip] 'NarrowPeak)
                  , [SCATACSeq S (File '[Gzip] 'Other)] )
               -> ReaderT config IO (Maybe (FilePath, FilePath, FilePath))
computePeakRAS _ (Nothing, _) = return Nothing
computePeakRAS prefix (peakFl, inputs) = do
    dir <- asks ((<> asDir prefix) . _scatacseq_output_dir) >>= getPath
    let output1 = dir <> "relative_accessibility_scores.tsv"
        output2 = dir <> "cell_specificity_score.tsv"
        output3 = dir <> "/cell_specificity_pvalue.tsv"
    liftIO $ do
        peaks <- fmap (map mkName) $ runResourceT $ runConduit $
            streamBedGzip (fromJust peakFl^.location) .| sinkList

        (names, cols) <- fmap unzip $ forM inputs $ \input -> do
            vec <- computeRAS (input^.replicates._2.files.location)
            return (input^.eid, V.convert vec)

        let ras = DF.fromMatrix peaks names $ Mat.fromColumns cols
            ss = computeSS ras
        DF.writeTable output1 (T.pack . show) ras
        DF.writeTable output2 (T.pack . show) ss
        cdf <- computeCDF ras
        DF.writeTable output3 (T.pack . show) $ DF.map (lookupP cdf) ss
        return $ Just (output1, output2, output3)
  where
    mkName :: BED3 -> T.Text
    mkName p = T.pack $ B.unpack (p^.chrom) <> ":" <> show (p^.chromStart) <>
        "-" <> show (p^.chromEnd)
    lookupP (vec, res, n) x | p == 0 = 1 / n
                            | otherwise = p
      where
        p = vec V.! i
        i = min (V.length vec - 1) $ truncate $ x / res 