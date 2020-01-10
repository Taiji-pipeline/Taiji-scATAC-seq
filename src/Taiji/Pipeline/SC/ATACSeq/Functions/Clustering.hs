{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE DataKinds #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Clustering
    ( filterMatrix
    , spectral
    , mkKNNGraph
    , clustering
        
    , plotClusters
    , extractTags
    , subMatrix
    , splitBedByCluster'
    ) where

import qualified Data.ByteString.Char8 as B
import Data.Binary (encodeFile, decodeFile)
import Bio.Utils.Functions (scale)
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
import Taiji.Pipeline.SC.ATACSeq.Functions.Feature (streamMatrices)
import qualified Taiji.Utils.DataFrame as DF
import Taiji.Utils.Plot
import Taiji.Utils
import Taiji.Utils.Plot.ECharts

filterMatrix :: (Elem 'Gzip tags ~ 'True, SCATACSeqConfig config)
             => FilePath
             -> SCATACSeq S (File tags 'Other)
             -> ReaderT config IO (SCATACSeq S (File '[] 'Tsv, File tags 'Other))
filterMatrix prefix input = do
    dir <- asks ((<> asDir prefix) . _scatacseq_output_dir) >>= getPath
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
         -> SCATACSeq S (a, File tags 'Other)
         -> ReaderT config IO (SCATACSeq S (a, File '[Gzip] 'Tsv))
spectral prefix seed input = do
    dir <- asks ((<> asDir prefix) . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_spectral.tsv.gz" dir
            (T.unpack $ input^.eid) (input^.replicates._1)
    input & replicates.traversed.files %%~ liftIO . ( \(rownames, fl) -> do
        shelly $ run_ "taiji-utils" $ ["reduce", T.pack $ fl^.location,
            T.pack output] ++ maybe []
            (\x -> ["--seed", T.pack $ show x]) seed
        return (rownames, location .~ output $ emptyFile)
        )

mkKNNGraph :: SCATACSeqConfig config
           => FilePath
           -> SCATACSeq S (File '[] 'Tsv, File '[Gzip] 'Tsv)
           -> ReaderT config IO (SCATACSeq S (File '[] 'Tsv, File '[] 'Other, File '[] Tsv))
mkKNNGraph prefix input = do
    dir <- asks ((<> asDir ("/" ++ prefix)) . _scatacseq_output_dir) >>= getPath
    let output_knn = printf "%s/%s_rep%d_knn.npz" dir (T.unpack $ input^.eid)
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
    tmp <- asks _scatacseq_tmp_dir
    dir <- asks ((<> asDir ("/" ++ prefix)) . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_clusters.bin" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traversed.files %%~ liftIO . ( \(idx, knn, umap) -> withTemp tmp $ \tmpFl -> do
        shelly $ run_ "taiji-utils" $
            [ "clust", T.pack $ knn^.location, T.pack tmpFl
            , "--res", T.pack $ show resolution
            , "--optimizer"
            , case optimizer of
                RBConfiguration -> "RB"
                CPM -> "CPM"
            ]
        let sourceCells = getZipSource $ (,,) <$>
                ZipSource (iterateC succ 0) <*>
                ZipSource ( sourceFile (idx^.location) .|
                    linesUnboundedAsciiC .|
                    mapC ((\[a,b] -> (a,b)) . B.split '\t') ) <*>
                ZipSource ( sourceFile (umap^.location) .|
                    linesUnboundedAsciiC .|
                    mapC (map readDouble . B.split '\t') )
        cells <- runResourceT $ runConduit $ sourceCells .| mapC f .| sinkVector
        clusters <- map (map readInt . B.split ',') . B.lines <$>
            B.readFile tmpFl
        let cellCluster = zipWith (\i -> CellCluster $ B.pack $ "C" ++ show i) [1::Int ..] $
                map (map (cells V.!)) clusters
        encodeFile output cellCluster
        return $ location .~ output $ emptyFile )
  where
    f (i, (bc, dep), [d1,d2]) = Cell i (d1,d2) bc $ readInt dep
    f _ = error "formatting error"

-- | Extract tags for clusters.
extractTags :: FilePath   -- ^ Directory to save the results
            -> Builder ()
extractTags prefix = do
    nodePar "Extract_Tags" 'splitBedByCluster $ return ()
    uNode "Merge_Tags_Prep" [| \input ->
        map (first head . unzip) $ groupBy ((==) `on` fst) $
        sortBy (comparing fst) $ concat $ input^..folded.replicates._2.files
        |]
    nodePar "Merge_Tags" [| mergeBedCluster prefix |] $ return ()
    path ["Extract_Tags", "Merge_Tags_Prep", "Merge_Tags"]

-- | Extract BEDs for each cluster.
splitBedByCluster :: SCATACSeqConfig config
                  => ( SCATACSeq S (File '[NameSorted, Gzip] 'Bed, a)
                     , SCATACSeq S (File '[] 'Other) ) -- ^ Cluster files
                  -> ReaderT config IO
                       (SCATACSeq S [(B.ByteString, File '[] 'Bed)])
splitBedByCluster (input, clFl) = do
    let idRep = asDir $ "/temp/" <> T.unpack (input^.eid) <>
            "_rep" <> show (input^.replicates._1)
        exp_id = B.pack $ T.unpack $ input^.eid
    dir <- asks _scatacseq_output_dir >>= getPath . (<> idRep)
    input & replicates.traverse.files %%~ ( \(bed,_) -> liftIO $ do
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
          -> [SCATACSeq S (File tags 'Other)]   -- ^ Matrices
          -> File tag' 'Other   -- Clusters
          -> ReaderT config IO [SCATACSeq S (File tags 'Other)]
subMatrix prefix mats clFl = do
    dir <- asks _scatacseq_output_dir >>= getPath . (<> (asDir prefix))
    liftIO $ do
        cls <- decodeFile $ clFl^.location
        mat <- mkSpMatrix id $ head mats ^. replicates._2.files.location
        let mkSink CellCluster{..} = filterC ((`S.member` ids) . fst) .|
                (sinkRows (S.size ids) (_num_col mat) id output >> return res)
              where
                ids = S.fromList $ map _cell_barcode _cluster_member
                output = dir <> B.unpack _cluster_name <> ".mat.gz"
                res = head mats & eid .~ T.pack (B.unpack _cluster_name)
                                & replicates._2.files.location .~ output
        runResourceT $ runConduit $ streamMatrices id mats .|
            sequenceSinks (map mkSink cls)


--------------------------------------------------------------------------------
-- Vizualization
--------------------------------------------------------------------------------

plotClusters :: FilePath
             -> (FilePath, SCATACSeq S (File '[] 'Other))
             -> IO ()
plotClusters dir (qc, input) = do
    stats <- readStats qc
    inputData <- decodeFile $ input^.replicates._2.files.location
    let output = printf "%s/%s_rep%d_cluster.html" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
        output2 = printf "%s/%s_metadata.tsv" dir (T.unpack $ input^.eid)
        (nms, num_cells) = unzip $ map (\(CellCluster nm cells) ->
            (T.pack $ B.unpack nm, fromIntegral $ length cells)) inputData
        plt = stackBar $ DF.mkDataFrame ["number of cells"] nms [num_cells]
        compos = composition inputData
    clusters <- sampleCells inputData
    savePlots output [] $ visualizeCluster clusters ++
        clusterComposition compos : tissueComposition compos : plt :
            clusterQC stats inputData
    outputMetaData output2 stats inputData
  where
    clusterQC stats cls =
        [ plotNumReads res
        , plotTE res
        , plotDoubletScore res
        , plotDupRate res ]
      where
        res = flip map cls $ \x ->
            (T.pack $ B.unpack $ _cluster_name x, map f $ _cluster_member x)
        f x = M.lookupDefault undefined (_cell_barcode x) statMap
        statMap = M.fromList $ map (\x -> (_barcode x, x)) stats

outputMetaData :: FilePath -> [Stat] -> [CellCluster] -> IO ()
outputMetaData output stat = B.writeFile output . B.unlines . (header:) . concatMap f
  where
    f CellCluster{..} = map g _cluster_member
      where
        g Cell{..} = B.intercalate "\t"
            [ _cell_barcode
            , cl
            , toShortest $ fst _cell_2d
            , toShortest $ snd _cell_2d
            , toShortest $ maybe 0 _te $ M.lookup _cell_barcode stat'
            , fromJust $ packDecimal $ maybe 0 _uniq_reads $ M.lookup _cell_barcode stat'
            , toShortest $ maybe 0 _doublet_score $ M.lookup _cell_barcode stat'
            , toShortest $ maybe 0 _mito_rate $ M.lookup _cell_barcode stat'
            , toShortest $ maybe 0 _dup_rate $ M.lookup _cell_barcode stat'
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