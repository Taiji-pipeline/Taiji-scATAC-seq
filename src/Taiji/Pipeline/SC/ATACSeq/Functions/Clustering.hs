{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE DeriveLift #-}
{-# LANGUAGE StandaloneDeriving #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Clustering
    ( plotClusters
    , spectralClust
    , extractTags
    , extractSubMatrix
    , subMatrix
    , getBedCluster
    , extractBedByBarcode 
    , ClustOpt(..)
    , toParams
    , defClustOpt
    , mkKNNGraph
    , clustering
    , combineClusters
    ) where

import qualified Data.ByteString.Char8 as B
import Data.Binary (encodeFile, decodeFile)
import Language.Haskell.TH.Syntax (Lift)
import qualified Data.Conduit.List as CL
import qualified Data.Text as T
import qualified Data.HashSet as S
import qualified Data.HashMap.Strict as M
import Data.Singletons.Prelude (Elem)
import Bio.Data.Bed
import Control.Arrow (first)
import qualified Data.Vector as V
import System.IO
import Data.List.Ordered (nubSort)
import Bio.Utils.Misc (readDouble, readInt)
import Shelly (shelly, run_, escaping)
import Control.Workflow
import Data.Conduit.Zlib (gzip)
   
import Taiji.Pipeline.SC.ATACSeq.Functions.QC
import Taiji.Pipeline.SC.ATACSeq.Functions.DR
import Taiji.Prelude
import Taiji.Pipeline.SC.ATACSeq.Types
import Taiji.Pipeline.SC.ATACSeq.Functions.Feature (streamMatrices)
import qualified Taiji.Utils.DataFrame as DF
import Taiji.Utils.Plot
import Taiji.Utils
import Taiji.Utils.Plot.ECharts

-- | Perform spectral clustering.
spectralClust :: FilePath   -- ^ Directory to save the results
              -> ClustOpt -> Builder ()
spectralClust prefix opt = do
    nodePar "Filter_Mat" [| filterMatrix prefix |] $ return ()
    nodePar "Reduce_Dims" [| spectral prefix Nothing |] $ return ()
    nodePar "Make_KNN" [| \x -> mkKNNGraph prefix $
        x & replicates.traverse.files %~ return
        |] $ return ()
    node "Cluster_Config" [| \xs -> if null xs
        then return []
        else do
            dir <- asks _scatacseq_output_dir >>= getPath . (<> (asDir prefix))
            let output = dir ++ "/parameters.txt"
            liftIO $ writeFile output $ unlines $ map
                (\x -> T.unpack (x^.eid) ++ "\t" ++ show (_resolution opt)) xs
            return $ zip (repeat output) xs
        |] $ return ()
    nodePar "Cluster" [| \(fl, x) -> do
        let f [a,b] = (T.pack a, read b)
        ps <- liftIO $ map (f . words) . lines <$> readFile fl
        clustering prefix opt{_resolution = fromJust $ lookup (x^.eid) ps} x 
        |] $ return ()
    path ["Filter_Mat", "Reduce_Dims", "Make_KNN", "Cluster_Config", "Cluster"]

-- | Clustering options
data ClustOpt = ClustOpt
    { _optimizer :: Optimizer 
    , _resolution :: Double
    } deriving (Lift)

deriving instance Lift Optimizer

defClustOpt :: ClustOpt
defClustOpt = ClustOpt RBConfiguration 1 

toParams :: ClustOpt -> [T.Text]
toParams ClustOpt{..} = 
    [ "--res"
    , T.pack $ show _resolution
    , "--optimizer"
    , case _optimizer of
        RBConfiguration -> "RB"
        CPM -> "CPM"
    ]

mkKNNGraph :: SCATACSeqConfig config
           => FilePath
           -> SCATACSeq S [(File '[] 'Tsv, File '[Gzip] 'Tsv)]
           -> ReaderT config IO (SCATACSeq S (File '[] 'Tsv, File '[] 'Other, File '[] Tsv))
mkKNNGraph prefix input = do
    dir <- asks ((<> asDir ("/" ++ prefix)) . _scatacseq_output_dir) >>= getPath
    let output_knn = printf "%s/%s_rep%d_knn.npz" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
        output_umap = printf "%s/%s_rep%d_umap.txt" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traversed.files %%~ liftIO . ( \fls -> do
        shelly $ run_ "taiji-utils" $
            [ "knn"
            , T.pack $ intercalate "," $ map (^.location) $ snd $ unzip fls
            , T.pack output_knn
            , "-k", "50"
            , "--embed", T.pack output_umap
            , "--thread", "4" ]
        return ( head $ fst $ unzip fls
               , location .~ output_knn $ emptyFile
               , location .~ output_umap $ emptyFile )
        )

clustering :: SCATACSeqConfig config
           => FilePath
           -> ClustOpt
           -> SCATACSeq S (File '[] 'Tsv, File '[] 'Other, File '[] Tsv)
           -> ReaderT config IO (SCATACSeq S (File '[] 'Other))
clustering prefix opt input = do
    tmp <- asks _scatacseq_tmp_dir
    dir <- asks ((<> asDir ("/" ++ prefix)) . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_clusters.bin" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traversed.files %%~ liftIO . ( \(idx, knn, umap) -> withTemp tmp $ \tmpFl -> do
        shelly $ run_ "taiji-utils" $
            ["clust", T.pack $ knn^.location, T.pack tmpFl] ++
            toParams opt

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
    f (i, (bc, dep), [d1,d2,d3,d4,d5]) = Cell i (d1,d2) (d3,d4,d5) bc $ readInt dep
    f _ = error "formatting error"

-- | Extract tags for clusters.
extractTags :: FilePath   -- ^ Directory to save the results
            -> Builder ()
extractTags prefix = do
    nodePar "Extract_Tags" 'getBedCluster $ return ()
    node "Merge_Tags_Prep" [| \input -> return $
        map (first head . unzip) $ groupBy ((==) `on` fst) $
        sortBy (comparing fst) $ concat $ input^..folded.replicates._2.files
        |] $ return ()
    nodePar "Merge_Tags" [| mergeBedCluster prefix |] $ return ()
    path ["Extract_Tags", "Merge_Tags_Prep", "Merge_Tags"]

-- | Extract BEDs for each cluster.
getBedCluster :: SCATACSeqConfig config
              => ( SCATACSeq S (File '[NameSorted, Gzip] 'Bed, a)
                 , SCATACSeq S (File '[] 'Other) ) -- ^ Cluster files
              -> ReaderT config IO
                    (SCATACSeq S [(B.ByteString, File '[] 'Bed)])
getBedCluster (input, clFl) = do
    let idRep = asDir $ "/temp/" <> T.unpack (input^.eid) <>
            "_rep" <> show (input^.replicates._1)
    dir <- asks _scatacseq_output_dir >>= getPath . (<> idRep)
    input & replicates.traverse.files %%~ ( \(bed,_) -> liftIO $ do
        clusters <- getClusterBarcodes (B.pack $ T.unpack $ input^.eid) clFl
        let outputs = M.fromList $
                map (\x -> (x, dir ++ "/" ++ B.unpack x ++ ".bed")) $
                nubSort $ M.elems clusters
        fileHandles <- mapM (flip openFile WriteMode) outputs
        runResourceT $ runConduit $ streamBedGzip (bed^.location) .|
            mapM_C ( \x -> case M.lookup (fromJust $ x^.name) clusters of
                Nothing -> return ()
                Just bc -> liftIO $
                    B.hPutStrLn (M.lookupDefault undefined bc fileHandles) $
                    toLine (x :: BED)
            )
        return $ M.toList $ fmap (\x -> location .~ x $ emptyFile) outputs )
  where
    getClusterBarcodes :: B.ByteString    -- ^ experiment id
                    -> SCATACSeq S (File '[] 'Other)   -- ^ Cluster files
                    -> IO (M.HashMap B.ByteString B.ByteString)  -- ^ Barcode to cluster map
    getClusterBarcodes exp_id cl = do
        clusters <- decodeFile $ cl^.replicates._2.files.location :: IO [CellCluster]
        return $ M.fromListWith (error "same barcode") $
            flip concatMap clusters $ \CellCluster{..} ->
                zip (mapMaybe getBarcode _cluster_member) $ repeat _cluster_name
      where
        getBarcode x | B.init i == exp_id = Just bc
                     | otherwise = Nothing
          where
            (i, bc) = B.breakEnd (=='+') $ _cell_barcode x

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

extractBedByBarcode :: (Elem 'Gzip tags ~ 'True, Elem 'NameSorted tags ~ 'True)
                    => [FilePath]   -- ^ Outputs
                    -> [[B.ByteString]]  -- ^ Barcodes
                    -> File tags 'Bed
                    -> IO [File '[] 'Bed]
extractBedByBarcode outputs bcs input = do
    fileHandles <- V.fromList <$> mapM (\x -> openFile x WriteMode) outputs
    let bcIdx = M.fromListWith (error "same barcode") $ concat $
            zipWith (\i x -> zip x $ repeat i) [0..] bcs
        f x = let idx = M.lookupDefault (error $ show nm) nm bcIdx
                  nm = fromJust ((head x :: BED) ^. name) 
              in liftIO $ B.hPutStrLn (fileHandles V.! idx) $ B.unlines $ map toLine x 
    runResourceT $ runConduit $ streamBedGzip (input^.location) .|
        CL.groupBy ((==) `on` (^.name)) .| mapM_C f
    return $ map (\x -> location .~ x $ emptyFile) outputs

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

-- | Extract submatrix
extractSubMatrix :: SCATACSeqConfig config
                 => FilePath   -- ^ dir
                 -> SCATACSeq S (File tags 'Other, File '[] 'Other)
                 -> ReaderT config IO [SCATACSeq S (File tags 'Other)]
extractSubMatrix prefix input = do
    dir <- asks _scatacseq_output_dir >>= getPath . (<> (asDir prefix))
    let (matFl, clFl) = input^.replicates._2.files
    liftIO $ do
        clusters <- decodeFile $ clFl^.location
        mat <- mkSpMatrix id $ matFl^.location
        let mkSink CellCluster{..} = filterC ((`S.member` ids) . fst) .|
                (yield header >> mapC (encodeRowWith id)) .| unlinesAsciiC .|
                gzip .| (sinkFile output >> return res)
              where
                res = input & eid .~ input^.eid <> "+" <> T.pack (B.unpack _cluster_name)
                            & replicates._2.files .~ (location .~ output $ matFl)
                output = dir <> T.unpack (input^.eid) <> "_" <> B.unpack _cluster_name <> ".mat.gz"
                header = B.pack $ printf "Sparse matrix: %d x %d" nCell
                    (_num_col mat)
                ids = S.fromList $ map _cell_barcode _cluster_member
                nCell = S.size ids
        runResourceT $ runConduit $ streamRows mat .|
            sequenceSinks (map mkSink clusters)

combineClusters :: SCATACSeqConfig config
                => FilePath
                -> [SCATACSeq S (File tags 'Other)]
                -> ReaderT config IO (SCATACSeq S (File tags 'Other))
combineClusters _ [a] = return a
combineClusters prefix inputs = do
    dir <- asks _scatacseq_output_dir >>= getPath . (<> (asDir prefix))
    let output = dir ++ "combined_clusters.bin"
    cls <- liftIO $ fmap concat $ forM inputs $ \x -> do
        let nm = B.pack $ T.unpack $ x^.eid 
        cs <- decodeFile $ x^.replicates._2.files.location :: IO [CellCluster]
        return $ flip map cs $ \c -> c{ _cluster_name= nm <> _cluster_name c }
    liftIO $ encodeFile output cls
    return $ head inputs & eid .~ "All"
                         & replicates._2.files.location .~ output
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
        (nms, num_cells) = unzip $ map (\(CellCluster nm cells) ->
            (T.pack $ B.unpack nm, fromIntegral $ length cells)) inputData
        plt = stackBar $ DF.mkDataFrame ["number of cells"] nms [num_cells]
        compos = composition inputData
    clusters <- sampleCells inputData
    savePlots output [] $ visualizeCluster clusters ++
        clusterComposition compos : tissueComposition compos : plt :
            clusterQC stats inputData
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