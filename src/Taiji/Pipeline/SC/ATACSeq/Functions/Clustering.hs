{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE DataKinds #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Clustering
    ( filterMatrix
    , spectral
    , batchCorrection
    , mkKNNGraph
    , clustering
    , clusterComposition
    , tissueComposition
    , composition
        
    , plotClusters
    , extractTags
    , subMatrix
    , splitBedByCluster
    , splitBedByCluster'
    ) where

import qualified Data.ByteString.Char8 as B
import Data.Binary (encodeFile, decodeFile)
import Bio.Utils.Functions (scale)
import Data.Conduit.Zlib (multiple, ungzip, gzip)
import qualified Data.Conduit.List as CL
import qualified Data.Text as T
import qualified Data.HashSet as S
import qualified Data.HashMap.Strict as M
import Data.Singletons.Prelude (Elem)
import Bio.Data.Bed
import Control.Arrow (first, (&&&))
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed.Mutable as UM
import qualified Data.Vector.Unboxed as U
import System.IO
import Data.List.Ordered (nubSort)
import Shelly (shelly, run_, escaping)
import Control.Workflow
   
import Taiji.Pipeline.SC.ATACSeq.Functions.QC
import Taiji.Prelude
import Taiji.Pipeline.SC.ATACSeq.Types
import qualified Taiji.Utils.DataFrame as DF
import Taiji.Utils
import Taiji.Utils.Plot (savePlots)
import Taiji.Utils.Plot.ECharts

filterMatrix :: (Elem 'Gzip tags ~ 'True, SCATACSeqConfig config)
             => FilePath
             -> SCATACSeq S (File tags 'Matrix)
             -> ReaderT config IO (SCATACSeq S (File '[] 'Tsv, File tags 'Matrix))
filterMatrix dir input = do
    let output = printf "%s/%s_rep%d_filt.mat.gz" dir
            (T.unpack $ input^.eid) (input^.replicates._1)
        rownames = printf "%s/%s_rep%d_rownames.txt" dir
            (T.unpack $ input^.eid) (input^.replicates._1)
    input & replicates.traversed.files %%~ ( \fl -> liftIO $ do
        sp <- mkSpMatrix readInt $ fl^.location
        runResourceT $ runConduit $
            streamRows sp .| mapC f .| unlinesAsciiC .| sinkFile rownames
        counts <- do
            v <- UM.replicate (_num_col sp) 0
            runResourceT $ runConduit $ streamRows sp .| concatMapC snd .|
                mapC fst .| mapM_C (UM.unsafeModify v (+1))
            U.unsafeFreeze v
        let (zeros, nonzeros) = U.partition ((==0) . snd) $
                U.zip (U.enumFromN 0 (U.length counts)) counts
            (i, v) = U.unzip nonzeros
            idx = U.toList $ fst $ U.unzip $ U.filter ((>1.65) . snd) $ U.zip i $ scale v
        filterCols output (idx ++ U.toList (fst $ U.unzip zeros)) $ fl^.location
        return ( location .~ rownames $ emptyFile
               , location .~ output $ emptyFile ) )
  where
    f (nm, xs) = nm <> "\t" <> fromJust (packDecimal $ foldl1' (+) $ map snd xs)

-- | Reduce dimensionality using spectral clustering
spectral :: (Elem 'Gzip tags ~ 'True, SCATACSeqConfig config)
         => FilePath  -- ^ directory
         -> Maybe Int  -- ^ seed
         -> SCATACSeq S (a, File tags 'Matrix)
         -> ReaderT config IO (SCATACSeq S (a, File '[Gzip] 'Tsv))
spectral dir seed input = do
    let output = printf "%s/%s_rep%d_spectral.tsv.gz" dir
            (T.unpack $ input^.eid) (input^.replicates._1)
    input & replicates.traversed.files %%~ liftIO . ( \(rownames, fl) -> do
        shelly $ run_ "taiji-utils" $ ["reduce", T.pack $ fl^.location,
            T.pack output] ++ maybe []
            (\x -> ["--seed", T.pack $ show x]) seed
        return (rownames, location .~ output $ emptyFile)
        )

batchCorrection :: SCATACSeqConfig config
                => String -> SCATACSeq S (File '[] 'Tsv, File '[Gzip] 'Tsv)
                -> ReaderT config IO (SCATACSeq S (File '[] 'Tsv, File '[Gzip] 'Tsv))
batchCorrection prefix input = do
    dir <- asks ((<> asDir prefix) . _scatacseq_output_dir) >>= getPath
    asks _scatacseq_batch_info >>= \case
        Nothing -> return input
        Just batchFl -> do
            let output = printf "%s/%s_rep%d_spectral_corrected.tsv.gz" dir
                    (T.unpack $ input^.eid) (input^.replicates._1)
            input & replicates.traversed.files %%~ liftIO . ( \(rownames, fl) -> do
                idToBatchMap <- M.fromListWith undefined <$> readBatchInfo batchFl
                let f x = let (i, r) = B.breakEnd (=='_') x
                          in case M.lookup (B.init i) idToBatchMap of
                              Nothing -> Nothing
                              Just (l, g) -> Just (l <> r, g)
                labels <- map (f . B.init . fst . B.breakEnd (=='+') . head . B.split '\t') . B.lines <$>
                    B.readFile (rownames^.location)
                if (all isNothing labels)
                    then return (rownames, fl)
                    else do
                        readData (fl^.location) >>= batchCorrect labels >>= writeData output
                        return (rownames, location .~ output $ emptyFile)
                )
  where
    readData fl = runResourceT $ runConduit $
        sourceFile fl .| multiple ungzip .| linesUnboundedAsciiC .|
        mapC (U.fromList . map readDouble . B.split '\t') .| sinkVector
    writeData output vec = runResourceT $ runConduit $ yieldMany vec .|
        mapC (B.intercalate "\t" . map toShortest . U.toList) .|
        unlinesAsciiC .| gzip .| sinkFile output


mkKNNGraph :: SCATACSeqConfig config
           => FilePath
           -> SCATACSeq S (File '[] 'Tsv, File '[Gzip] 'Tsv)
           -> ReaderT config IO (SCATACSeq S (File '[] 'Tsv, File '[] 'Other, File '[] Tsv))
mkKNNGraph dir input = do
    let output_knn = printf "%s/%s_rep%d_knn.txt.gz" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
        output_umap = printf "%s/%s_rep%d_umap.txt" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traversed.files %%~ liftIO . ( \(idxFl, matFl) -> do
        shelly $ run_ "taiji-utils" $
            [ "knn"
            , T.pack $ matFl^.location
            , T.pack output_knn
            , "-k", "50"
            , "--embed", T.pack output_umap
            , "--thread", "4" ]
        return ( idxFl
               , location .~ output_knn $ emptyFile
               , location .~ output_umap $ emptyFile )
        )

clustering :: SCATACSeqConfig config
           => FilePath
           -> Double
           -> Optimizer
           -> SCATACSeq S (File '[] 'Tsv, File '[] 'Other, File '[] Tsv)
           -> ReaderT config IO (SCATACSeq S (File '[] 'Other))
clustering prefix resolution optimizer input = do
    dir <- asks ((<> asDir ("/" ++ prefix)) . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_clusters.bin" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traversed.files %%~ liftIO . ( \(idx, knn, umap) -> do
        let sourceCells = getZipSource $ (,,) <$>
                ZipSource (iterateC succ 0) <*>
                ZipSource (sourceFile (idx^.location) .| linesUnboundedAsciiC .| mapC g) <*>
                ZipSource ( sourceFile (umap^.location) .|
                    linesUnboundedAsciiC .|
                    mapC (map readDouble . B.split '\t') )
        cells <- runResourceT $ runConduit $ sourceCells .| mapC f .| sinkVector
        let clusterSizeCutoff = max minCells $ min 50 $ V.length cells `div` 1000
        clusters <- fmap (sortBy (flip (comparing length)) . filter ((>=clusterSizeCutoff) . length)) $
            readKNNGraph (knn^.location) >>= leiden resolution optimizer
        let cellCluster = zipWith (\i x -> CellCluster (B.pack $ "C" ++ show i) x Nothing) [1::Int ..] $
                map (map (cells V.!)) clusters
        encodeFile output cellCluster
        return $ location .~ output $ emptyFile )
  where
    minCells = 10
    f (i, (bc, dep), [d1,d2]) = Cell i (d1,d2) bc $ readInt dep
    f _ = error "formatting error"
    g x = case B.split '\t' x of
        [a, b] -> (a, b)
        [a] -> (a, "0")
        _ -> error "formatting error"

-- | Extract tags for clusters.
extractTags :: FilePath   -- ^ Directory to save the results
            -> Builder ()
extractTags prefix = do
    nodePar "Extract_Tags" 'splitBedByCluster $ return ()
    uNode "Merge_Tags_Prep" [| \input -> return $
        map (first head . unzip) $ groupBy ((==) `on` fst) $
        sortBy (comparing fst) $ concat $ input^..folded.replicates._2.files
        |]
    nodePar "Merge_Tags" [| mergeBedCluster prefix |] $ return ()
    path ["Extract_Tags", "Merge_Tags_Prep", "Merge_Tags"]

-- | Extract BEDs for each cluster.
splitBedByCluster :: SCATACSeqConfig config
                  => ( SCATACSeq S (File '[NameSorted, Gzip] 'Bed)
                     , SCATACSeq S (File '[] 'Other) ) -- ^ Cluster files
                  -> ReaderT config IO
                       (SCATACSeq S [(B.ByteString, File '[] 'Bed)])
splitBedByCluster (input, clFl) = do
    let idRep = asDir $ "/temp/" <> T.unpack (input^.eid) <>
            "_rep" <> show (input^.replicates._1)
        exp_id = B.pack $
            T.unpack (input^.eid) <> "_" <> show (input^.replicates._1)
    dir <- asks _scatacseq_output_dir >>= getPath . (<> idRep)
    input & replicates.traverse.files %%~ ( \bed -> liftIO $ do
        let filterBarcode cl =
                let x = filter (B.isPrefixOf exp_id . _cell_barcode) $
                        _cluster_member cl
                in cl{_cluster_member=x}
            addId = name %~ fmap (\x -> exp_id <> "+" <> x)
        clusters <- fmap (map filterBarcode) $ decodeFile $
            clFl^.replicates._2.files.location
        runResourceT $ runConduit $ streamBedGzip (bed^.location) .|
            mapC addId .| splitBedByCluster' dir clusters
        )

-- | Extract BEDs for each cluster.
splitBedByCluster' :: MonadIO m
                  => FilePath   -- ^ Output directory
                  -> [CellCluster]
                  -> ConduitT BED o m [(B.ByteString, File '[] 'Bed)]
splitBedByCluster' dir clusters = do
    let (nm, bcs) = unzip $
            map (_cluster_name &&& map _cell_barcode . _cluster_member) clusters 
        bcIdxMap = M.fromListWith (error "same barcode") $ concat $
            zipWith (\x -> zip x . repeat) bcs [0..]
        outputs = map (\x -> dir <> "/" <> B.unpack x <> ".bed") nm
    fileHandles <- liftIO $ V.fromList <$>
        mapM (\x -> openFile x WriteMode) outputs
    CL.groupBy ((==) `on` (^.name)) .| mapM_C ( \beds ->
        maybe (return ())
            (\x -> liftIO $ B.hPutStr (fileHandles V.! x) $ B.unlines $ map toLine beds) $
            M.lookup (fromJust $ head beds ^. name) bcIdxMap
            )
    liftIO $ mapM_ hClose fileHandles
    return $ zip nm $ map (\x -> location .~ x $ emptyFile) outputs
{-# INLINE splitBedByCluster' #-}

-- | Extract BEDs for each cluster.
mergeBedCluster :: SCATACSeqConfig config
                => FilePath  -- ^ prefix
                -> (B.ByteString, [File '[] 'Bed])
                -> ReaderT config IO (B.ByteString, File '[Gzip] 'Bed)
mergeBedCluster prefix (cName, fls) = do
    dir <- asks _scatacseq_output_dir >>= getPath . (<> asDir prefix)
    let output = dir <> B.unpack cName <> ".bed.gz"
    liftIO $ do
        shelly $ escaping False $ do
            run_ "cat" $ map (\x -> T.pack $ x^.location) fls ++
                ["|", "gzip", "-c", ">", T.pack output]
            run_ "rm" $ map (\x -> T.pack $ x^.location) fls
        return (cName, location .~ output $ emptyFile)

-- | Extract cluster submatrix
subMatrix :: SCATACSeqConfig config
          => FilePath   -- ^ Dir
          -> [SCATACSeq S (File tags 'Matrix)]   -- ^ Matrices
          -> File tag' 'Other   -- Clusters
          -> ReaderT config IO [SCATACSeq S (File tags 'Matrix)]
subMatrix prefix inputs clFl = do
    dir <- asks _scatacseq_output_dir >>= getPath . (<> (asDir prefix))
    liftIO $ do
        cls <- decodeFile $ clFl^.location
        mat <- mkSpMatrix id $ head inputs ^. replicates._2.files.location
        mergedMat <- fmap concatMatrix $ forM inputs $ \input -> do
            let f x = B.concat [ B.pack $ T.unpack $ input^.eid, "_"
                    , B.pack $ show (input^.replicates._1), "+", x ]
            fmap (mapRows (first f)) $ mkSpMatrix id $ input^.replicates._2.files.location
        let mkSink CellCluster{..} = filterC ((`S.member` ids) . fst) .|
                (sinkRows (S.size ids) (_num_col mat) id output >> return res)
              where
                ids = S.fromList $ map _cell_barcode _cluster_member
                output = dir <> B.unpack _cluster_name <> ".mat.gz"
                res = head inputs & eid .~ T.pack (B.unpack _cluster_name)
                                  & replicates._2.files.location .~ output
        runResourceT $ runConduit $ streamRows mergedMat .| sequenceSinks (map mkSink cls)

--------------------------------------------------------------------------------
-- Vizualization
--------------------------------------------------------------------------------

plotClusters :: SCATACSeqConfig config
             => ( [(Double, (Int, Double, Double))]
                , [(Double, ([FilePath], FilePath))]
                , Maybe (SCATACSeq S (File '[] 'Tsv, File '[] 'Other, File '[] Tsv)) )
             -> ReaderT config IO (Maybe (SCATACSeq S (File '[] 'Other)))
plotClusters (_, _, Nothing) = return Nothing
plotClusters (params, clFl, Just knn) = do
    dir <- asks _scatacseq_output_dir >>= getPath . (<> "/Cluster/")
    resAuto <- liftIO $ optimalParam (dir <> "/Clustering_parameters.html") params
    res <- fromMaybe resAuto <$> asks _scatacseq_cluster_resolution 
    let clOutput = printf "%s/%s_rep%d_clusters.bin" dir (T.unpack $ knn^.eid)
            (knn^.replicates._1)
    fmap Just $ knn & replicates.traversed.files %%~ liftIO . ( \(idx, _, umap) -> do
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
        encodeFile clOutput cellCluster
        let clViz = printf "%s/%s_rep%d_cluster.html" dir (T.unpack $ knn^.eid)
                (knn^.replicates._1)
            (nms, num_cells) = unzip $ map (\(CellCluster nm cs _) ->
                (T.pack $ B.unpack nm, fromIntegral $ length cs)) cellCluster
            plt = stackBar $ DF.mkDataFrame ["number of cells"] nms [num_cells]
            compos = composition cellCluster
        clusters' <- sampleCells cellCluster
        savePlots clViz [] $ visualizeCluster clusters' ++
            figRepro : clusterComposition compos : tissueComposition compos : [plt]
        return $ location .~ clOutput $ emptyFile )
  where
    f (i, (bc, dep), [d1,d2]) = Cell i (d1,d2) bc $ readInt dep
    f _ = error "formatting error"
    g x = case B.split '\t' x of
        [a, b] -> (a, b)
        [a] -> (a, "0")
        _ -> error "formatting error"

-- | Compute the normalized tissue composition for each cluster.
tissueComposition :: DF.DataFrame Int -> EChart
tissueComposition = stackBar . DF.map round' . DF.mapCols normalize .
    DF.transpose . DF.mapCols normalize . DF.map fromIntegral
  where
    round' x = fromIntegral (round $ x * 1000 :: Int) / 1000
    normalize :: V.Vector Double -> V.Vector Double
    normalize xs | V.all (==0) xs = xs
                 | otherwise = V.map (/V.sum xs) xs

-- | Compute the cluster composition for each tissue.
clusterComposition :: DF.DataFrame Int -> EChart
clusterComposition = stackBar . DF.map round' . DF.mapCols normalize .
    DF.map fromIntegral
  where
    round' x = fromIntegral (round $ x * 1000 :: Int) / 1000
    normalize :: V.Vector Double -> V.Vector Double
    normalize xs | V.all (==0) xs = xs
                 | otherwise = V.map (/V.sum xs) xs

-- | Compute the cluster x tissue composition table.
composition :: [CellCluster] -> DF.DataFrame Int
composition clusters = DF.mkDataFrame rownames colnames $
    flip map rows $ \x -> map (\i -> M.lookupDefault 0 i x) colnames
  where
    (rownames, rows) = unzip $ flip map clusters $ \CellCluster{..} ->
        ( T.pack $ B.unpack _cluster_name
        , M.fromListWith (+) $ map (\x -> (getName x, 1)) _cluster_member )
    colnames = nubSort $ concatMap M.keys rows
    getName Cell{..} =
        let prefix = fst $ B.breakEnd (=='+') _cell_barcode
        in if B.null prefix then "" else T.pack $ B.unpack $ B.init prefix