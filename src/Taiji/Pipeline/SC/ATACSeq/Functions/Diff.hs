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
    , rpkmPeak
    , rpkmDiffPeak
    , plotDiffGene
    ) where

import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U
import qualified Data.ByteString.Char8 as B
import qualified Data.Text as T
import qualified Data.Text.IO as T
import Language.Javascript.JMacro
import Bio.Data.Bed
import Bio.Utils.Functions (scale)
import Control.Arrow (second)
import Data.Binary (decodeFile)
import Bio.Data.Bed.Utils (rpkmBed)
import Data.Conduit.Internal (zipSources)
import System.Random.MWC (create)
import System.Random.MWC.Distributions (uniformShuffle)
import Data.List.Ordered (nubSort)
import qualified Data.HashSet as S
import qualified Data.HashMap.Strict as M
import qualified Data.Matrix.Unboxed as MU
import Data.Conduit.Zlib (gzip)
import Shelly hiding (FilePath)

import Taiji.Prelude
import Taiji.Pipeline.SC.ATACSeq.Types
import Taiji.Pipeline.SC.ATACSeq.Functions.Utils
import qualified Taiji.Utils.DataFrame as DF
import Taiji.Utils.Plot
import qualified Taiji.Utils.Plot.ECharts as E

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
         -> Bool
         -> SCATACSeq S ([[B.ByteString]], File '[Gzip] 'Other)
         -> ReaderT config IO (SCATACSeq S (File '[Gzip] 'Other))
mkRefMat prefix addName input = do
    dir <- asks ((<> asDir prefix) . _scatacseq_output_dir) >>= getPath
    let output = dir <> "/" <> T.unpack (input^.eid) <> "_ref_mat.gz"
    input & replicates.traversed.files %%~ ( \(bcs, fl) -> liftIO $ do
        mat <- mkSpMatrix id $ fl^.location
        let header = B.pack $ printf "Sparse matrix: %d x %d" (S.size bc) (_num_col mat)
            bc = S.fromList $ concat bcs
            filterFn x | addName = (B.pack (T.unpack $ input^.eid) <> "+" <> x) `S.member` bc
                       | otherwise = x `S.member` bc
        runResourceT $ runConduit $ streamRows mat .| filterC (filterFn . fst) .|
            (yield header >> mapC (encodeRowWith id)) .| unlinesAsciiC .|
            gzip .| sinkFile output
        return $ location .~ output $ emptyFile )

diffPeaks :: SCATACSeqConfig config
          => ( SCATACSeq S (File '[Gzip] 'Bed, File tags 'Other)
             , File '[Gzip] 'Other ) -- ^ Ref matrix
          -> ReaderT config IO (SCATACSeq S (File '[Gzip] 'NarrowPeak))
diffPeaks (input, ref) = do
    dir <- asks ((<> "/Diff/Peak/") . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d.np.gz" dir
            (T.unpack $ input^.eid) (input^.replicates._1)
    input & replicates.traversed.files %%~ liftIO . ( \(peakFl, fl) -> do
        stats <- fmap (M.fromList . map (\(a,b,c,d) -> (a, (b,c,d))) . filter (\x -> x^._4 < 0.01)) $
            diffAnalysis (fl^.location) (ref^.location) Nothing
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


-- | Compute RPKM for each peak
rpkmPeak :: SCATACSeqConfig config
         => ( (B.ByteString, File '[Gzip] 'Bed)
            , Maybe (File '[Gzip] 'NarrowPeak) )
         -> ReaderT config IO (B.ByteString, FilePath)
rpkmPeak ((nm, tags), Just peakFl) = do
    dir <- asks _scatacseq_output_dir >>= getPath . (<> "/Feature/Peak/")
    let output = dir <> "peak_rpkm_" <> B.unpack nm <> ".txt"
    liftIO $ do
        peaks <- runResourceT $ runConduit $
            streamBedGzip (peakFl^.location) .| sinkList :: IO [BED3]
        vec <- runResourceT $ runConduit $ streamBedGzip (tags^.location) .| rpkmBed peaks 
        B.writeFile output $ B.unlines $ map toShortest $ U.toList vec
        return (nm, output)
rpkmPeak _ = undefined

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
          -> Maybe FilePath
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
             => [SCATACSeq S (File '[] 'Tsv)]
             -> ReaderT config IO ()
plotDiffGene inputs = do
    dir <- figDir
    (cls, fdrs) <- liftIO $ fmap unzip $ forM inputs $ \input -> do
        fdr <- readFDR $ input^.replicates._2.files.location
        return (input^.eid, fdr)
    markers <- asks _scatacseq_marker_gene_list >>= \case
        Nothing -> return []
        Just fl -> liftIO $ readMarkers fl
    let output = dir <> "/diff_gene.html"
        genes = nubSort $ concatMap M.keys fdrs
        df1 = DF.mkDataFrame cls (map fst markers) $ flip map fdrs $ \fdr ->
            U.toList $ scale' $ U.fromList $ map snd $ enrichment markers fdr
        df2 = DF.mkDataFrame cls genes $ flip map fdrs $ \fdr ->
            map (\g -> M.lookupDefault 0 g fdr) genes
    liftIO $ savePlots output [] [mkHeatmap df1, mkHeatmap df2]
  where
    mkHeatmap df = p <> E.toolbox <> E.option
        [jmacroE| { visualMap: { inRange: {color: `viridis`} } } |]
      where
        p = E.heatmap $ DF.reorderRows (DF.orderByCluster id) $
            DF.reorderColumns (DF.orderByCluster id) df
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
             -> Maybe FilePath  -- ^ selected idx
             -> IO [(Int, Double, Double, Double)]
diffAnalysis fg bk idx = withTemp Nothing $ \tmp -> do
    shelly $ run_ "sc_utils" $ [ "diff", T.pack tmp, "--fg", T.pack fg
        , "--bg", T.pack bk ] ++ fromMaybe [] (fmap (\x -> ["--index", T.pack x]) idx)
    map (f . B.words) . B.lines <$> B.readFile tmp
  where
    f [a,b,c,d] = (readInt a, readDouble b, readDouble c, readDouble d)
    f _ = error "wrong format!"