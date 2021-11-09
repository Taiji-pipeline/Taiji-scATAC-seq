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
    , vizCluster
    , plotSubclusters
    ) where

import qualified Data.ByteString.Char8 as B
import qualified Data.Map.Strict as Map
import Data.Conduit.Zlib (multiple, ungzip, gzip)
import qualified Data.Text as T
import Data.Conduit.Internal (zipSources)
import qualified Data.HashMap.Strict as M
import qualified Data.Vector as V
import qualified Data.HashSet as S
import Data.Binary (decodeFile, encodeFile)
import qualified Data.Vector.Unboxed as U
import Shelly (shelly, run_)
   
import Taiji.Pipeline.SC.ATACSeq.Functions.QC
import Taiji.Pipeline.SC.ATACSeq.Functions.Clustering
import qualified Taiji.Utils.DataFrame as DF
import Taiji.Prelude
import Taiji.Pipeline.SC.ATACSeq.Types
import Taiji.Utils
import Taiji.Utils.Plot (savePlots)
import Taiji.Utils.Plot.ECharts


-- | Reduce dimensionality using spectral clustering
subSpectral :: SCATACSeqConfig config
            => SCATACSeq S (File '[] 'Tsv, File '[Gzip] 'Matrix)
            -> ReaderT config IO (SCATACSeq S (File '[] 'Tsv, File '[Gzip] 'Tsv))
subSpectral input = do
    dir <- asks ((<> "/Subcluster/Spectra/") . _scatacseq_output_dir) >>= getPath
    tmpdir <- asks _scatacseq_tmp_dir
    let output = printf "%s/%s_spectra.tsv.gz" dir (T.unpack $ input^.eid)
    withRunInIO $ \runInIO -> withTempDir tmpdir $ \tmp -> runInIO $
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
                            let tmp' = tmp <> "/tmp.tsv.gz"
                            shelly $ run_ "taiji-utils" $ ["reduce", T.pack $ fl^.location,
                                T.pack tmp', "--seed", "23948"]
                            readData tmp' >>= batchCorrect labels >>= writeData output
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
                 , [ SCATACSeq S ( File '[RowName, Gzip] 'Tsv
                               , File '[ColumnName, Gzip] 'Tsv
                               , File '[Gzip] 'Matrix) ] )
              -> ReaderT config IO (SCATACSeq S (File '[] 'Tsv, File '[Gzip] 'Matrix))
subsetFeatMat (clFl, matFls) = do
    let prefix = asDir $ "/Subcluster/Matrix/" <> T.unpack (clFl^.eid)
    tmp <- asks _scatacseq_tmp_dir
    dir <- asks ((<> prefix) . _scatacseq_output_dir) >>= getPath
    let clId = B.pack $ T.unpack $ clFl^.eid
    [cluster] <- liftIO $ fmap (filter ((==clId) . _cluster_name)) $ decodeFile $
        clFl^.replicates._2.files.location
    let barcodes = S.fromList $ map _cell_barcode $ _cluster_member cluster
    withRunInIO $ \runInIO -> withTempDir tmp $ \tmpdir -> runInIO $
        ( clFl & replicates.traversed.files %%~ ( \_ -> liftIO $ do
            mat <- fmap concatMatrix $ forM matFls $ \x -> mkSpMatrix id $ x^.replicates._2.files._3.location
            runResourceT $ runConduit $ streamRows mat .|
                filterC ((`S.member` barcodes) . fst) .|
                sinkRows (S.size barcodes) (_num_col mat) id (tmpdir <> "/mat.gz")
            return $ location .~ (tmpdir <> "/mat.gz") $ emptyFile
            ) ) >>= filterMatrix dir

combineClusters :: SCATACSeqConfig config
                => ( Maybe (SCATACSeq S (File '[] 'Other))
                   , [SCATACSeq S (File '[] 'Other)] )
                -> ReaderT config IO (Maybe (SCATACSeq S (File '[] 'Other)))
combineClusters (cl, []) = return cl
combineClusters (Just cl, subCl) = do
    exclude <- S.fromList <$> asks _scatacseq_cluster_exclude
    dir <- asks ((<> "/Cluster/") . _scatacseq_output_dir) >>= getPath
    let output = dir <> "final_clusters.bin"
        clIds = S.fromList $ map (B.pack . T.unpack . (^.eid)) subCl
    liftIO $ do
        cls <- forM subCl $ \e -> do
            let i = B.pack $ T.unpack $ e^.eid
            c <- decodeFile $ e^.replicates._2.files.location
            return $ if length c == 1
                then [(head c){_cluster_name = i}]
                else map (rename (B.pack $ T.unpack $ e^.eid)) $
                    filter ((>=0.5) . fromMaybe 1 . _cluster_reproducibility) c
        cs <- fmap (filter (not . (`S.member` clIds) . _cluster_name)) $
            decodeFile $ cl^.replicates._2.files.location
        encodeFile output $
            filter (not . (`S.member` exclude) . T.pack . B.unpack . _cluster_name) $
            concat cls <> cs
    return $ Just $ replicates._2.files.location .~ output $ cl
  where
    rename prefix c = c{_cluster_name = prefix <> "." <> B.tail (_cluster_name c)}
combineClusters _ = return Nothing

vizCluster :: SCATACSeqConfig config
           => ( Maybe (SCATACSeq S (File '[] 'Other))
              , Maybe (SCATACSeq S (File '[] 'Tsv, File '[Gzip] 'Tsv))
              , FilePath )
           -> ReaderT config IO ()
vizCluster (Just clFl, Just coord, qc) = do
    tmp <- asks _scatacseq_tmp_dir
    fig <- figDir
    liftIO $ do
        cellCluster <- withTempDir tmp $ \tmpdir -> do
            cls <- decodeFile $ clFl^.replicates._2.files.location
            let s1 = sourceFile (coord^.replicates._2.files._1.location) .|
                    linesUnboundedAsciiC .| mapC (head . B.split '\t')
                s2 = sourceFile (coord^.replicates._2.files._2.location) .|
                    multiple ungzip .| linesUnboundedAsciiC
            runResourceT $ runConduit $ zipSources s1 s2 .| umap tmpdir cls

        stats <- readStats qc
        let clViz = fig <> "/final_cluster.html" 
            clMeta = fig <> "/cell_metadata.tsv"
            (nms, num_cells) = unzip $ map (\(CellCluster nm cs _) ->
                (T.pack $ B.unpack nm, fromIntegral $ length cs)) cellCluster
            plt = stackBar $ DF.mkDataFrame ["number of cells"] nms [num_cells]
            compos = composition cellCluster
        clusters' <- sampleCells cellCluster
        savePlots clViz [] $ visualizeCluster clusters' ++
            clusterComposition compos : tissueComposition compos : plt : []
        outputMetaData clMeta stats cellCluster
  where
    clusterQC stats cls = statToJson res
      where
        res = flip map cls $ \x ->
            (T.pack $ B.unpack $ _cluster_name x, map h $ _cluster_member x)
        h x = M.lookupDefault
            (error $ "barcode not found: " <> show (_cell_barcode x)) (_cell_barcode x) statMap
        statMap = M.fromList $ map (\x -> (_barcode x, x)) stats
vizCluster _ = return ()

umap :: FilePath -> [CellCluster]
     -> ConduitT (B.ByteString, B.ByteString) Void (ResourceT IO) [CellCluster]
umap dir clusters = do
    (bcs, _, _) <- concatMapC f .| sink
    coord <- liftIO $ do
        shelly $ run_ "taiji-utils" ["viz", T.pack outData, T.pack out, "-t", T.pack outCl]
        map (map readDouble . B.split '\t') . B.lines <$> B.readFile out
    let coord' = M.fromList $ zip bcs coord
        updateCoord c = let [x, y] = M.lookupDefault undefined (_cell_barcode c) coord'
                        in c{_cell_2d = (x,y)}
    return $ flip map clusters $ \cl -> cl{_cluster_member = map updateCoord $ _cluster_member cl}
  where
    sink = getZipSink $ (,,) <$>
        ZipSink (mapC (^._1) .| sinkList) <*>
        ZipSink (mapC (^._2) .| unlinesAsciiC .| sinkFile outCl) <*>
        ZipSink (mapC (^._3) .| unlinesAsciiC .| sinkFile outData)
    outCl = dir <> "/cl.txt"
    outData = dir <> "/data.txt"
    out = dir <> "/umap.txt"
    f (bc, x) = case M.lookup bc bcMap of
        Nothing -> Nothing
        Just c -> Just (bc, c, x)
    bcMap = M.fromList $ flip concatMap clusters $ \cl ->
        zip (map _cell_barcode $ _cluster_member cl) $ repeat $ getClusterName cl
    --getClusterName x = fst $ B.break (=='.') $ _cluster_name x
    getClusterName x = _cluster_name x
{-# INLINE umap #-}

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
    knn & replicates.traversed.files %%~ liftIO . ( \(idx, _, u) -> do
        let sourceCells = getZipSource $ (,,) <$>
                ZipSource (iterateC succ 0) <*>
                ZipSource (sourceFile (idx^.location) .| linesUnboundedAsciiC .| mapC g) <*>
                ZipSource ( sourceFile (u^.location) .|
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

outputMetaData :: FilePath -> [Stat] -> [CellCluster] -> IO ()
outputMetaData output stat = B.writeFile output . B.unlines . (header:) . concatMap f
  where
    f CellCluster{..} = map g _cluster_member
      where
        g Cell{..} = B.intercalate "\t" $ map (fromMaybe "NA")
            [ Just _cell_barcode
            , Just cl
            , Just $ toShortest $ fst _cell_2d
            , Just $ toShortest $ snd _cell_2d
            , fmap (toShortest . _te) $ M.lookup _cell_barcode stat'
            , fmap (fromJust . packDecimal . _uniq_reads) $ M.lookup _cell_barcode stat'
            , join $ fmap (fmap toShortest . _doublet_score) $ M.lookup _cell_barcode stat'
            , join $ fmap (fmap toShortest . _mito_rate) $ M.lookup _cell_barcode stat'
            , join $ fmap (fmap toShortest . _dup_rate) $ M.lookup _cell_barcode stat'
            ]
        cl = B.tail _cluster_name
    header = B.intercalate "\t"
        [ "Sample+Barcode"
        , "Cluster"
        , "UMAP1"
        , "UMAP2"
        , "TSSe"
        , "Num_fragments"
        , "Doublet_score"
        , "Fraction_Mito"
        , "Fraction_duplication"
        ]
    stat' = M.fromList $ map (\x -> (_barcode x, x)) stat
{-# INLINE outputMetaData #-}

