{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE QuasiQuotes #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Diff
    ( sampleCells
    , mkRefMat
    , diffPeaks
    , diffGenes
    , rpkmDiffPeak
    , plotDiffGene
    , specificPeaks
    ) where

import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U
import qualified Data.ByteString.Char8 as B
import qualified Data.Text as T
import qualified Data.Text.IO as T
import Language.Javascript.JMacro
import Bio.Data.Bed
import Bio.Utils.Functions (scale, filterFDR, binarySearch)
import Control.Arrow (second)
import Data.Binary (decodeFile)
import Data.Conduit.Internal (zipSources)
import System.Random.MWC (create)
import System.Random.MWC.Distributions (uniformShuffle)
import Data.List.Ordered (nubSort)
import qualified Data.HashSet as S
import qualified Data.HashMap.Strict as M
import qualified Data.Matrix.Unboxed as MU
import qualified Data.Matrix as Mat
import Shelly hiding (FilePath)
import qualified Data.Vector.Algorithms.Intro as I

import Taiji.Prelude
import Taiji.Pipeline.SC.ATACSeq.Types
import Taiji.Pipeline.SC.ATACSeq.Functions.Utils
import qualified Taiji.Utils.DataFrame as DF
import Taiji.Utils.Plot
import qualified Taiji.Utils.Plot.ECharts as E
import Taiji.Pipeline.SC.ATACSeq.Functions.Feature (streamMatrices)

sampleCells :: Int   -- ^ Number of cells
            -> SCATACSeq S (File tags 'Other)   -- ^ clusters
            -> IO (SCATACSeq S [[B.ByteString]])
sampleCells n input = input & replicates.traversed.files %%~ ( \fl -> do
    cls <- decodeFile $ fl^.location
    g <- create
    forM cls $ \cl -> do
        let v = V.fromList $ map _cell_barcode $ _cluster_member cl
        V.toList . V.take n <$> uniformShuffle v g
    )

mkRefMat :: SCATACSeqConfig config
         => FilePath   -- ^ Prefix
         -> [B.ByteString]   -- ^ barcodes
         -> [SCATACSeq S (File '[Gzip] 'Other)]    -- ^ matrices
         -> ReaderT config IO (Maybe (File '[Gzip] 'Other))
mkRefMat _ _ [] = return Nothing
mkRefMat prefix bcs mats = do
    dir <- asks ((<> asDir prefix) . _scatacseq_output_dir) >>= getPath
    let output = dir <> "/ref_mat.gz"
    liftIO $ do
        mat <- mkSpMatrix id $ head mats ^. replicates._2.files.location
        runResourceT $ runConduit $ streamMatrices id mats .|
            filterC ((`S.member` bcSet) . fst) .| sinkRows (S.size bcSet) (_num_col mat) id output
        return $ Just $ location .~ output $ emptyFile
  where
    bcSet = S.fromList bcs

diffPeaks :: SCATACSeqConfig config
          => FilePath 
          -> ( File '[Gzip] 'NarrowPeak
             , SCATACSeq S (File tags 'Other, File '[] 'NarrowPeak)   -- ^ Cluster peaks
             , File '[Gzip] 'Other ) -- ^ Ref matrix
          -> ReaderT config IO (SCATACSeq S (File '[Gzip] 'NarrowPeak))
diffPeaks prefix (peakFl, input, ref) = do
    dir <- asks ((<> asDir prefix) . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d.np.gz" dir
            (T.unpack $ input^.eid) (input^.replicates._1)
    input & replicates.traversed.files %%~ liftIO . ( \(fl, clPeak) -> do
        masterPeakList <- runResourceT $ runConduit $
            streamBedGzip (peakFl^.location) .| sinkList :: IO [BED3]
        peaks <- bedToTree const . map (\x -> (x, ())) <$>
            (readBed $ clPeak^.location :: IO [BED3])
        let idx = map fst $ filter (\(_, x) -> peaks `isIntersected` x) $
                zip [0..] masterPeakList
        stats <- fmap (M.fromList . map (\(a,b,c,d) -> (a, (b,c,d))) . filter (\x -> x^._4 < 0.01)) $
            diffAnalysis (fl^.location) (ref^.location) $ Just idx
        let f :: (Int, BED3) -> Maybe NarrowPeak
            f (i, peak) = case M.lookup i stats of
                Nothing -> Nothing
                Just (fold, pval, fdr) -> 
                    let logP | pval == 0 = 200
                             | otherwise = negate $ logBase 10 pval 
                        logFDR | fdr == 0 = 200
                               | otherwise = negate $ logBase 10 fdr
                    in Just $ npSignal .~ fold $
                        npPvalue .~ Just logP $
                        npQvalue .~ Just logFDR $ convert peak
        runResourceT $ runConduit $
            zipSources (iterateC succ 0) (streamBedGzip $ peakFl^.location) .|
            concatMapC f .| sinkFileBedGzip output
        return $ location .~ output $ emptyFile )

rpkmDiffPeak :: SCATACSeqConfig config
             => ( [SCATACSeq S (File '[Gzip] 'NarrowPeak)]   -- ^ Diff peaks
                , Maybe (File '[Gzip] 'NarrowPeak)  -- ^ All peaks
                , [(B.ByteString, FilePath)] )  -- ^ RPKM
             -> ReaderT config IO FilePath
rpkmDiffPeak (diffs, Just peak, rpkms) = do
    dir <- asks ((<> "/Diff/Peak/") . _scatacseq_output_dir) >>= getPath
    let output = dir <> "diff_peak_RPKM.tsv"
    liftIO $ do
        mat <- readRPKMs $ map snd rpkms
        peakList <- runResourceT $ runConduit $
            streamBedGzip (peak^.location) .| sinkList :: IO [BED3]
        res <- forM diffs $ \input -> do
            peaks <- fmap (filter f) $ runResourceT $ runConduit $
                streamBedGzip (input^.replicates._2.files.location) .| sinkList
            let idx = getPeakIndex peaks peakList
                submat = map (B.intercalate "\t" . map toShortest . U.toList) $
                    MU.toRows $ getSubMatrix idx mat
            return $ "#" <> B.pack (T.unpack $ input^.eid) : submat
        B.writeFile output $ B.unlines $ header : concat res
        return output
  where
    f x = x^.npSignal > 2
    header = B.intercalate "\t" $ fst $ unzip rpkms
    getSubMatrix idx mat = MU.fromRows $ map (mat `MU.takeRow`) idx
    readRPKMs fls = MU.fromColumns <$> mapM readData fls
      where
        readData fl = U.fromList . map readDouble . B.lines <$> B.readFile fl
rpkmDiffPeak _ = undefined

-- | Get the indices of query peaks in the reference peak list.
getPeakIndex :: (BEDLike b1, BEDLike b2)
             => [b1]   -- ^ Differential peaks
             -> [b2]   -- ^ Reference peaks
             -> [Int]  -- ^ Indices
getPeakIndex query ref = flip map query $ \q -> M.lookupDefault undefined
    (q^.chrom, q^.chromStart, q^.chromEnd) peakIdxMap
  where
    peakIdxMap = M.fromList $ zipWith
        (\x i -> ((x^.chrom, x^.chromStart, x^.chromEnd), i)) ref [0..]

diffGenes :: SCATACSeqConfig config
          => FilePath
          -> Maybe [Int]
          -> ( FilePath   -- ^ gene list
             , SCATACSeq S (File tags 'Other) 
             , File '[Gzip] 'Other ) -- ^ Ref matrix
          -> ReaderT config IO (SCATACSeq S (File '[] 'Tsv))
diffGenes prefix idx (nameFl, input, ref) = do
    dir <- asks ((<> asDir prefix) . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d.tsv" dir
            (T.unpack $ input^.eid) (input^.replicates._1)
    input & replicates.traversed.files %%~ liftIO . ( \fl -> do
        stats <- fmap (M.fromList . map (\(a,b,c,d) -> (a, (b,c,d))) . filter (\x -> x^._4 < 0.01)) $
            diffAnalysis (fl^.location) (ref^.location) idx
        let f (i, nm) = case M.lookup i stats of
                Nothing -> Nothing
                Just (fold, pval, fdr) ->
                    let logP | pval == 0 = 200
                             | otherwise = negate $ logBase 10 pval 
                        logFDR | fdr == 0 = 200
                               | otherwise = negate $ logBase 10 fdr
                    in Just $ B.intercalate "\t" [nm, toShortest fold,
                            toShortest logP, toShortest logFDR]
        runResourceT $ runConduit $ zipSources (iterateC succ 0)
            (sourceFile nameFl .| linesUnboundedAsciiC .| mapC (head . B.words)) .|
            concatMapC f .| unlinesAsciiC .| sinkFile output
        return $ location .~ output $ emptyFile )

plotDiffGene :: SCATACSeqConfig config
             => FilePath
             -> [SCATACSeq S (File '[] 'Tsv)]
             -> ReaderT config IO ()
plotDiffGene filename inputs = do
    dir <- figDir
    (cls, fdrs) <- liftIO $ fmap unzip $ forM inputs $ \input -> do
        fdr <- readFDR $ input^.replicates._2.files.location
        return (input^.eid, fdr)
    markers <- asks _scatacseq_marker_gene_list >>= \case
        Nothing -> return []
        Just fl -> liftIO $ readMarkers fl
    let output = dir <> filename
        genes = nubSort $ concatMap M.keys fdrs
        df1 = DF.mkDataFrame cls (map fst markers) $ flip map fdrs $ \fdr ->
            U.toList $ scale' $ U.fromList $ map snd $ enrichment markers fdr
        df2 = DF.mkDataFrame cls genes $ flip map fdrs $ \fdr ->
            map (\g -> logBase 2 $ M.lookupDefault 1 g fdr) genes
    liftIO $ do
        savePlots output [] [mkHeatmap df1]
        DF.writeTable (output <> ".tsv") (T.pack . show) df2
  where
    mkHeatmap df = p <> E.toolbox <> E.option
        [jmacroE| { visualMap: { inRange: {color: `viridis`} } } |]
      where
        p = E.heatmap $ DF.orderDataFrame id df
    scale' xs | U.all (==0) xs = xs
              | otherwise = scale xs
    enrichment markers fdr = map (second f) markers
      where
        f xs = foldl' (+) 0 (map (\k -> M.lookupDefault 0 k fdr) xs) /
            fromIntegral (length xs)
    readFDR fl = M.fromList . map snd . filter ((>2) . fst) .
        map (f . T.splitOn "\t") . T.lines <$> T.readFile fl
      where
        f (a:fld:_:q:_) = (read $ T.unpack fld :: Double, (a, read $ T.unpack q))
        f _ = undefined
    readMarkers :: FilePath -> IO [(T.Text, [T.Text])]
    readMarkers fl = M.toList . M.fromListWith (++) . map (f . T.splitOn "\t") .
        T.lines <$> T.readFile fl
      where
        f (a:b:_) = (b, [T.toUpper a])
        f _ = undefined

diffAnalysis :: FilePath  -- ^ foreground
             -> FilePath  -- ^ background
             -> Maybe [Int] -- ^ selected idx
             -> IO [(Int, Double, Double, Double)]
diffAnalysis fg bk idx = withTempDir Nothing $ \dir -> do
    let resFl = dir ++ "/res.txt"
        idxFl = dir ++ "/idx.txt"
        makeIdx x = do
            B.writeFile idxFl $ B.unlines $ map (B.pack . show) x
            return ["--index", T.pack idxFl]
    useIdx <- maybe (return []) makeIdx idx
    shelly $ run_ "sc_utils" $
        [ "diff", T.pack resFl, "--fg", T.pack fg
        , "--bg", T.pack bk ] ++ useIdx
    map (f . B.words) . B.lines <$> B.readFile resFl
  where
    f [a,b,c,d] = (readInt a, readDouble b, readDouble c, readDouble d)
    f _ = error "wrong format!"

-- | Get cell-specific peaks
specificPeaks :: SCATACSeqConfig config
             => FilePath
             -> ReaderT config IO [(T.Text, File '[Gzip] 'NarrowPeak)]
specificPeaks fl = do
    dir <- asks ((<> "/Diff/Peak/") . _scatacseq_output_dir) >>= getPath
    liftIO $ do
        df <- DF.readTable fl
        nullModel <- mkBackgroundModel $ DF._dataframe_data df
        let k = truncate $ n * fdr
            n = fromIntegral $ V.length nullModel
            bk = V.take k $ V.modify (\x -> I.partialSort x k) nullModel
            diff = Mat.toColumns $ Mat.fromRows $ map (snd . entropy) $
                Mat.toRows $ DF._dataframe_data df
            peaks = V.fromList $ map mkPeak $ DF.rowNames df
        forM (zip (DF.colNames df) diff) $ \(nm, xs) -> do
            let f i x | x > V.last bk = Nothing
                      | otherwise =
                          let p = fromIntegral (binarySearch bk x) / n
                          in Just $ npSignal .~ x $ npPvalue .~ Just p $ convert $ peaks V.! i
                output = dir <> T.unpack nm <> "_diff_peaks.np.gz"
                candidates = V.toList $ V.imapMaybe f xs
            runResourceT $ runConduit $ yieldMany candidates .| sinkFileBedGzip output
            return (nm, location .~ output $ emptyFile)
  where
    mkPeak :: T.Text -> BED3
    mkPeak x = let [chr, x'] = T.splitOn ":" x
                   [s,e] = T.splitOn "-" x'
               in asBed (B.pack $ T.unpack chr) (read $ T.unpack s) (read $ T.unpack e)
    fdr = 0.025

mkBackgroundModel :: Mat.Matrix Double -> IO (V.Vector Double)
mkBackgroundModel = fmap (V.concat . map (snd . entropy) . Mat.toRows) . shuffleMatrix
  where
    shuffleMatrix :: Mat.Matrix Double -> IO (Mat.Matrix Double)
    shuffleMatrix mat = fmap (Mat.fromVector (Mat.dim mat)) $
        create >>= uniformShuffle (Mat.flatten mat)

entropy :: V.Vector Double -> (Double, V.Vector Double)
entropy xs = (e, V.map (\x -> e - logBase 2 x) xs')
    where
    e = negate $ V.sum $ V.map f xs'
    xs' = V.map (/s) xs
    f p | p == 0 = 0
        | otherwise = p * logBase 2 p
    s = V.sum xs