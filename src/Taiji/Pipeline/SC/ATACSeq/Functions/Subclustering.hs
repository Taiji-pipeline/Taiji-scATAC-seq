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
import Data.Conduit.Zlib (multiple, ungzip, gzip)
import qualified Data.Text as T
import qualified Data.HashMap.Strict as M
import qualified Data.Vector as V
import qualified Data.HashSet as S
import Data.Binary (decodeFile, encodeFile)
import qualified Data.Vector.Unboxed as U
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

subsetFeatMat :: SCATACSeqConfig config
              => ( SCATACSeq S (File '[] 'Other)  -- ^ cluster
                 , SCATACSeq S (File '[Gzip] 'Other) )
              -> ReaderT config IO (SCATACSeq S (File '[] 'Tsv, File '[Gzip] 'Other))
subsetFeatMat (clFl, matFl) = do
    let prefix = asDir $ "/Subcluster/Matrix/" <> T.unpack (clFl^.eid)
    tmp <- asks _scatacseq_tmp_dir
    dir <- asks ((<> prefix) . _scatacseq_output_dir) >>= getPath
    let clId = B.pack $ T.unpack $ clFl^.eid
    [cluster] <- liftIO $ fmap (filter ((==clId) . _cluster_name)) $ decodeFile $
        clFl^.replicates._2.files.location
    let barcodes = S.fromList $ map _cell_barcode $ _cluster_member cluster
    withRunInIO $ \runInIO -> withTempDir tmp $ \tmpdir -> runInIO $
        ( clFl & replicates.traversed.files %%~ ( \_ -> liftIO $ do
            mat <- mkSpMatrix id $ matFl^.replicates._2.files.location
            runResourceT $ runConduit $ streamRows mat .|
                filterC ((`S.member` barcodes) . fst) .|
                sinkRows' (_num_col mat) id (tmpdir <> "/mat.gz")
            return $ location .~ (tmpdir <> "/mat.gz") $ matFl^.replicates._2.files
            ) ) >>= filterMatrix dir

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