{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE DataKinds #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Subclustering
    ( subSpectral
    , subsetFeatMat
    , combineClusters
    , plotSubclusters
    ) where

import qualified Data.ByteString.Char8 as B
import qualified Data.Map.Strict as Map
import Bio.Data.Bed hiding (NarrowPeak(..))
import Data.Conduit.Zlib (multiple, ungzip, gzip)
import qualified Data.Text as T
import qualified Data.HashMap.Strict as M
import qualified Data.Vector as V
import qualified Data.HashSet as S
import Data.Binary (decodeFile, encodeFile)
import qualified Data.Vector.Unboxed as U
import Data.List.Ordered (nubSort)
import Shelly (shelly, run_)
   
import Taiji.Pipeline.SC.ATACSeq.Functions.Clustering
import qualified Taiji.Utils.DataFrame as DF
import Taiji.Prelude
import Taiji.Pipeline.SC.ATACSeq.Types
import Taiji.Utils
import Taiji.Utils.Plot (savePlots)
import Taiji.Utils.Plot.ECharts


-- | Reduce dimensionality using spectral clustering
subSpectral :: SCATACSeqConfig config
            => SCATACSeq S (File '[] 'Tsv, File '[Gzip] 'Other)
            -> ReaderT config IO (SCATACSeq S (File '[] 'Tsv, File '[Gzip] 'Tsv))
subSpectral input = do
    dir <- asks ((<> "/Subcluster/Spectra/") . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_spectra.tsv.gz" dir (T.unpack $ input^.eid)
    withRunInIO $ \runInIO -> withTempDir (Just "./") $ \tmpdir -> runInIO $
        input & replicates.traversed.files %%~ ( \(rownames, fl) -> do
            asks _scatacseq_batch_info >>= \case
                Nothing -> liftIO $ shelly $ run_ "taiji-utils" $ ["reduce", T.pack $ fl^.location,
                    T.pack output, "--seed", "23948"]
                Just batchFl -> liftIO $ do
                    idToBatchMap <- M.fromListWith undefined <$> readBatchInfo batchFl
                    let f x = let (i, r) = B.breakEnd (=='_') x
                            in case M.lookup (B.init i) idToBatchMap of
                                    Nothing -> Nothing
                                    Just (l, g) -> Just (l <> r, g)
                    labels <- map (f . B.init . fst . B.breakEnd (=='+') . head . B.split '\t') . B.lines <$>
                        B.readFile (rownames^.location)
                    if (all isNothing labels)
                        then shelly $ run_ "taiji-utils" $ ["reduce", T.pack $ fl^.location,
                            T.pack output, "--seed", "23948"]
                        else do
                            let tmp = tmpdir <> "/tmp.tsv.gz"
                            shelly $ run_ "taiji-utils" $ ["reduce", T.pack $ fl^.location,
                                T.pack tmp, "--seed", "23948"]
                            readData tmp >>= batchCorrect labels >>= writeData output
            return ( rownames, location .~ output $ emptyFile )
            )
  where
    readData fl = runResourceT $ runConduit $
        sourceFile fl .| multiple ungzip .| linesUnboundedAsciiC .|
        mapC (U.fromList . map readDouble . B.split '\t') .| sinkVector
    writeData output vec = runResourceT $ runConduit $ yieldMany vec .|
        mapC (B.intercalate "\t" . map toShortest . U.toList) .|
        unlinesAsciiC .| gzip .| sinkFile output

{-
-- | subsetting of previous spectra
subSpectral :: SCATACSeqConfig config
            => SCATACSeq S (File '[] 'Tsv, File '[Gzip] 'Tsv, File '[] 'Other)
            -> ReaderT config IO (SCATACSeq S (File '[] 'Tsv, File '[Gzip] 'Tsv))
subSpectral input = do
    dir <- asks ((<> "/Subcluster/Spectra/") . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_spectra.tsv.gz" dir (T.unpack $ input^.eid)
        outNames = printf "%s/%s_rownames.tsv" dir (T.unpack $ input^.eid)
    input & replicates.traversed.files %%~ ( \(rownames, spec, clFl) -> do
        let clId = B.pack $ T.unpack $ input^.eid
        [cluster] <- liftIO $ fmap (filter ((==clId) . _cluster_name)) $ decodeFile $ clFl^.location
        let barcodes = S.fromList $ map _cell_barcode $ _cluster_member cluster
        _ <- runResourceT $ runConduit $
            zipSources (sourceFile (rownames^.location) .| linesUnboundedAsciiC)
                (sourceFile (spec^.location) .| multiple ungzip .| linesUnboundedAsciiC) .|
                filterC ((`S.member` barcodes) . head . B.split '\t' . fst) .|
                zipSinks (mapC fst .| unlinesAsciiC .| sinkFile outNames)
                    (mapC snd .| unlinesAsciiC .| gzip .| sinkFile output)
        return ( location .~ outNames $ emptyFile, location .~ output $ emptyFile )
        )
-}

subsetFeatMat :: SCATACSeqConfig config
              => ( SCATACSeq S (File '[] 'Other)  -- ^ cluster
                 --, [SCATACSeq S (File '[] 'Other)]  -- ^ original clusters
                 --, [SCATACSeq S [(B.ByteString, File '[Gzip] 'NarrowPeak)]]  -- ^ peaks
                 --, File '[Gzip] 'NarrowPeak
                 , SCATACSeq S (File '[Gzip] 'Other) )
              -> ReaderT config IO (SCATACSeq S (File '[] 'Tsv, File '[Gzip] 'Other))
subsetFeatMat (clFl, {-oriCl, peakFl, peakList,-} matFl) = do
    --dir <- asks ((<> "/Subcluster/Matrix/") . _scatacseq_output_dir) >>= getPath
    let clId = B.pack $ T.unpack $ clFl^.eid
    [cluster] <- liftIO $ fmap (filter ((==clId) . _cluster_name)) $ decodeFile $
        clFl^.replicates._2.files.location
    --peaks <- liftIO $ findOriCluster cluster oriCl >>= getPeaks peakFl
    --liftIO $ runResourceT $ runConduit $ yieldMany peaks .| sinkFileBedGzip (dir <> T.unpack (clFl^.eid) <> "_feat.bed.gz")
    --r <- liftIO $ runResourceT $ runConduit $ streamBedGzip (peakList^.location) .|
    --    intersectBedWith ((\_ x -> x) :: BED3 -> [BED3] -> [BED3]) peaks .| sinkList
    --let idx = map snd $ filter (not . null . fst) $ zip r [0..] :: [Int]
    let barcodes = S.fromList $ map _cell_barcode $ _cluster_member cluster
    withRunInIO $ \runInIO -> withTempDir (Just "./") $ \tmpdir -> runInIO $
        ( clFl & replicates.traversed.files %%~ ( \_ -> liftIO $ do
            mat <- {-fmap (selectCols idx) $-} mkSpMatrix id $ matFl^.replicates._2.files.location
            runResourceT $ runConduit $ streamRows mat .|
                filterC ((`S.member` barcodes) . fst) .|
                sinkRows' (_num_col mat) id (tmpdir <> "/mat.gz")
            return $ location .~ (tmpdir <> "/mat.gz") $ matFl^.replicates._2.files
            ) ) >>= filterMatrix ("/Subcluster/Matrix/" <> T.unpack (clFl^.eid))

getPeaks :: [SCATACSeq S [(B.ByteString, File '[Gzip] 'NarrowPeak)]]
         -> [(T.Text, [B.ByteString])]
         -> IO [BED3]
getPeaks peakFls cls = do
    ps <- runResourceT $ runConduit $ mapM_ streamBedGzip fls .| sinkList
    runConduit $ mergeBed ps .| sinkList
  where
    fls = map (\x -> M.lookupDefault undefined x peakFls') cls'
    peakFls' = M.fromList $ flip concatMap peakFls $ \e -> 
        let i = e^.eid <> "_" <> T.pack (show $ e^.replicates._1)
        in flip map (e^.replicates._2.files) $ \(c, fl) -> ((i, c), fl^.location)
    cls' = nubSort $ flip concatMap cls $ \(i, xs) -> zip (repeat i) xs

findOriCluster :: CellCluster -> [SCATACSeq S (File '[] 'Other)] -> IO [(T.Text, [B.ByteString])]
findOriCluster cl ori_cls = fmap catMaybes $ forM ori_cls $ \e -> do
    let i = e^.eid <> "_" <> T.pack (show $ e^.replicates._1)
    if i `S.member` exps
        then do
            cls <- decodeFile $ e^.replicates._2.files.location
            let r = flip mapMaybe cls $ \c -> if any (`S.member` barcodes) $ map (addPrefix i . _cell_barcode) $ _cluster_member c
                    then Just $ _cluster_name c
                    else Nothing
            return $ Just (i, r)
        else return Nothing
  where
    exps = S.fromList $ map fst $ filter ((>0.01) . snd) $ normalize $ M.toList $
        M.fromListWith (+) $ zip (map (T.init . fst . T.breakOnEnd "+" .
            T.pack . B.unpack . _cell_barcode) $ _cluster_member cl) $ repeat (1 :: Double)
    barcodes = S.fromList $ map _cell_barcode $ _cluster_member cl
    addPrefix x y = B.pack (T.unpack x) <> "+" <> y
    normalize xs = let s = foldl1' (+) $ map snd xs in map (\(a,b) -> (a, b / s)) xs

combineClusters :: SCATACSeqConfig config
                => [SCATACSeq S (File '[] 'Other)]
                -> ReaderT config IO (Maybe (SCATACSeq S (File '[] 'Other)))
combineClusters [] = return Nothing
combineClusters input = do
    dir <- asks ((<> "/Subcluster/") . _scatacseq_output_dir) >>= getPath
    let output = dir <> "allClusters.bin"
    liftIO $ do
        cls <- forM input $ \e -> do
            let i = B.pack $ T.unpack $ e^.eid
            c <- decodeFile $ e^.replicates._2.files.location
            return $ if length c == 1
                then [(head c){_cluster_name = i}]
                else map (rename (B.pack $ T.unpack $ e^.eid)) $
                    filter ((>=0.5) . fromMaybe 1 . _cluster_reproducibility) c
        encodeFile output $ concat cls
    return $ Just $ eid .~ "Merged" $ replicates._2.files.location .~ output $ head input 
  where
    rename prefix cl = cl{_cluster_name = prefix <> "." <> B.tail (_cluster_name cl)}

plotSubclusters :: SCATACSeqConfig config
             => ( [(Double, (Int, Double, Double))]
                , [(Double, ([FilePath], FilePath))]
                , SCATACSeq S (File '[] 'Tsv, File '[] 'Other, File '[] Tsv) )
             -> ReaderT config IO (SCATACSeq S (File '[] 'Other))
plotSubclusters (params, clFl, knn) = do
    dir <- asks _scatacseq_output_dir >>= getPath . (<> "/Subcluster/")
    let outputHtml = printf "%s/%s_cluster.html" dir (T.unpack $ knn^.eid)
        outputParam = printf "%s/%s_parameters.html" dir (T.unpack $ knn^.eid)
        outputCl = printf "%s/%s_cluster.bin" dir (T.unpack $ knn^.eid)
    resAuto <- liftIO $ optimalParam outputParam params
    res <- Map.findWithDefault resAuto (knn^.eid) .
        fromMaybe Map.empty <$> asks _scatacseq_subcluster_resolution 
    knn & replicates.traversed.files %%~ liftIO . ( \(idx, _, umap) -> do
        let sourceCells = getZipSource $ (,,) <$>
                ZipSource (iterateC succ 0) <*>
                ZipSource (sourceFile (idx^.location) .| linesUnboundedAsciiC .| mapC g) <*>
                ZipSource ( sourceFile (umap^.location) .|
                    linesUnboundedAsciiC .|
                    mapC (map readDouble . B.split '\t') )
            (perturbed, cl) = fromJust $ lookup res clFl
        clusters <- decodeFile cl
        perturbedCls <- mapM decodeFile perturbed
        let (clNames, repro) = unzip $ zipWith (\i x -> (T.pack $ "C" <> show i, x)) [1::Int ..] $
                computeReproducibility clusters perturbedCls
            figRepro = line $ DF.mkDataFrame ["reproducibility"] clNames [repro]
        cells <- runResourceT $ runConduit $ sourceCells .| mapC f .| sinkVector
        let cellCluster = zipWith3 (\i x y -> CellCluster (B.pack $ "C" ++ show i) x $ Just y)
                [1::Int ..] (map (map (cells V.!)) clusters) repro
        encodeFile outputCl cellCluster

        let (nms, num_cells) = unzip $ map (\(CellCluster nm cs _) ->
                (T.pack $ B.unpack nm, fromIntegral $ length cs)) cellCluster
            plt = stackBar $ DF.mkDataFrame ["number of cells"] nms [num_cells]
            compos = composition cellCluster
        clusters' <- sampleCells cellCluster
        savePlots outputHtml [] $ visualizeCluster clusters' ++
            [figRepro, clusterComposition compos, tissueComposition compos, plt]
        return $ location .~ outputCl $ emptyFile )
  where
    f (i, (bc, dep), [d1,d2]) = Cell i (d1,d2) bc $ readInt dep
    f _ = error "formatting error"
    g x = case B.split '\t' x of
        [a, b] -> (a, b)
        [a] -> (a, "0")
        _ -> error "formatting error"