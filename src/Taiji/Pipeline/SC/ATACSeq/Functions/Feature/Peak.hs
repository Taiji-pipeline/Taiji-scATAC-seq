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
import Bio.Data.Bed hiding (NarrowPeak)
import qualified Bio.Data.Bed as Bed
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
import Taiji.Utils

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
          -> ReaderT config IO (B.ByteString, File '[Gzip] 'NarrowPeak)
findPeaks prefix opts (cName, bedFl) = do
    tmpdir <- asks _scatacseq_tmp_dir
    dir <- asks _scatacseq_output_dir >>= getPath . (<> asDir prefix)
    let output = dir ++ "/" ++ B.unpack cName ++ ".narrowPeak.gz" 
    asks _scatacseq_blacklist >>= \case
        Nothing -> liftIO $ do
            r <- callPeaks output bedFl Nothing opts
            return (cName, r)
        Just blacklist -> liftIO $ withTemp tmpdir $ \tmp -> do
            _ <- callPeaks tmp bedFl Nothing opts
            blackRegions <- readBed blacklist :: IO [BED3]
            let bedTree = bedToTree const $ map (\x -> (x, ())) blackRegions
            runResourceT $ runConduit $
                (streamBedGzip tmp :: ConduitT () Bed.NarrowPeak (ResourceT IO) ()) .|
                filterC (not . isIntersected bedTree) .| sinkFileBedGzip output
            return (cName, location .~ output $ emptyFile)

-- | Merge peaks
mergePeaks :: SCATACSeqConfig config
           => FilePath
           -> [(B.ByteString, File '[Gzip] 'NarrowPeak)]
           -> ReaderT config IO (Maybe (File '[Gzip] 'NarrowPeak))
mergePeaks _ [] = return Nothing
mergePeaks prefix input = do
    tmpdir <- asks _scatacseq_tmp_dir
    dir <- asks _scatacseq_output_dir >>= getPath . (<> asDir prefix)
    let output = dir <> "/merged.narrowPeak.gz" 
    liftIO $ withTemp tmpdir $ \tmp1 -> withTemp tmpdir $ \tmp2 -> do
        runResourceT $ runConduit $ mapM_ (streamBedGzip . (^._2.location)) input .|
            mapC resize .| sinkFileBed tmp1
        shelly $ escaping False $ bashPipeFail bash_ "cat" $
            [T.pack tmp1, "|", "sort", "-k1,1", "-k2,2n", "-k3,3n", ">", T.pack tmp2]
        let source = streamBed tmp2 .|
                mergeSortedBedWith iterativeMerge .| concatC
        runResourceT $ runConduit $ zipSources (iterateC succ (0 :: Int)) source .|
            mapC (\(i, p) -> name .~ Just ("p" <> B.pack (show i)) $ p) .|
            sinkFileBedGzip output
    return $ Just $ location .~ output $ emptyFile
  where
    iterativeMerge [] = []
    iterativeMerge peaks = bestPeak : iterativeMerge rest
      where
        rest = filter (\x -> sizeOverlapped x bestPeak == 0) peaks
        bestPeak = maximumBy (comparing (fromJust . (^.npPvalue))) peaks
    resize pk = chromStart .~ max 0 (summit - 250) $
        chromEnd .~ summit + 250 $
        npPeak .~ Just 250 $ pk
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
    getValue :: [Bed.NarrowPeak]  -- ^ query
            -> [BED3]  -- ^ Peak list
            -> [(Double, Double)]
    getValue query peakList = runIdentity $ runConduit $
        yieldMany peakList .| intersectBedWith f query .| sinkList
      where
        f _ [] = (0, 0)
        f _ [x] = (x^.npSignal, quantile (chiSquared 1) $ 10**(negate $ fromJust $ x^.npPvalue))
        f _ _ = undefined
getPeakEnrichment _ = undefined

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
            ss = computeSS $ DF.map (logBase 2 . (+1)) ras
        DF.writeTable output1 (T.pack . show) ras
        DF.writeTable output2 (T.pack . show) ss
        cdf <- computeCDF $ DF.map (logBase 2 . (+1)) ras
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