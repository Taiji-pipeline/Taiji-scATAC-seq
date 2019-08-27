{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE DeriveLift #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Clustering
    ( plotClusters
    , spectralClust
    , extractTags
    , extractSubMatrix
    , subMatrix
    , getBedCluster
    , extractBedByBarcode 
    , Embedding(..)
    , Normalization(..)
    , ClustOpt(..)
    , defClustOpt
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
import System.Random.MWC.Distributions
import System.Random.MWC
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
import Taiji.Pipeline.SC.ATACSeq.Functions.Utils
import Taiji.Pipeline.SC.ATACSeq.Functions.Feature (streamMatrices)
import qualified Taiji.Utils.DataFrame as DF
import Taiji.Utils.Plot
import Taiji.Utils.Plot.ECharts

-- | Perform spectral clustering.
spectralClust :: FilePath   -- ^ Directory to save the results
              -> ClustOpt -> Builder ()
spectralClust prefix opt = do
    nodePar "Filter_Mat" [| filterMatrix prefix |] $ return ()
    nodePar "Reduce_Dims" [| spectral prefix Nothing |] $ return ()
    node "Cluster_Config" [| \xs -> do
        dir <- asks _scatacseq_output_dir >>= getPath . (<> (asDir prefix))
        r <- case _resolution opt of
            Nothing -> return Nothing
            Just r -> liftIO $ do
                let output = dir ++ "/parameters.txt"
                writeFile output $ unlines $ map
                    (\x -> T.unpack (x^.eid) ++ "\t" ++ show r) xs
                return $ Just output
        return $ zip (repeat r) xs
        |] $ return ()
    nodePar "Cluster" [| \(para, x) -> case para of
        Nothing -> clustering prefix opt $ x & replicates.traverse.files %~ return
        Just fl -> do
            let f [a,b] = (T.pack a, read b)
            ps <- liftIO $ map (f . words) . lines <$> readFile fl
            clustering prefix opt{_resolution = lookup (x^.eid) ps} $
                x & replicates.traverse.files %~ return
        |] $ return ()
    path ["Filter_Mat", "Reduce_Dims", "Cluster_Config", "Cluster"]

-- | Embedding method
data Embedding = UMAP
               | TSNE
               | NoEmbedding
               deriving (Lift)

-- | Normalization method
data Normalization = Drop1st
                   | UnitBall
                   | None
                   deriving (Lift)

-- | Clustering options
data ClustOpt = ClustOpt
    { _normalization :: Normalization
    , _embedding_method :: Embedding
    , _dim :: Maybe Int  -- ^ How many dimensions to be used
    , _neighbors :: Int
    , _resolution :: Maybe Double
    } deriving (Lift)

defClustOpt :: ClustOpt
defClustOpt = ClustOpt None UMAP Nothing 50 Nothing

toParams :: ClustOpt -> [T.Text]
toParams ClustOpt{..} = embed ++ normalize ++ dim ++ res ++
    ["-k", T.pack $ show _neighbors]
  where
    embed = case _embedding_method of
        UMAP -> ["--embed-method", "umap"]
        TSNE -> ["--embed-method", "tsne"]
        NoEmbedding -> ["--embed-method", "none"]
    normalize = case _normalization of
        Drop1st -> ["--discard"]
        UnitBall -> ["--scale"]
        None -> []
    dim = case _dim of
        Nothing -> []
        Just d -> ["--dim", T.pack $ show d]
    res = case _resolution of
        Nothing -> []
        Just r -> ["--res", T.pack $ show r]
 
clustering :: SCATACSeqConfig config
           => FilePath
           -> ClustOpt
           -> SCATACSeq S [(File '[] 'Tsv, File '[Gzip] 'Tsv)]
           -> ReaderT config IO (SCATACSeq S (File '[] 'Other))
clustering prefix opt input = do
    tmp <- asks _scatacseq_temp_dir
    dir <- asks ((<> asDir ("/" ++ prefix)) . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_clusters.bin" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traversed.files %%~ liftIO . ( \fl -> do
        clustering' opt tmp fl >>= encodeFile output
        return $ location .~ output $ emptyFile )

clustering' :: ClustOpt
            -> Maybe FilePath      -- ^ temp dir
            -> [(File '[] 'Tsv, File '[Gzip] 'Tsv)]   -- ^ lsa input matrix
            -> IO [CellCluster]
clustering' opt dir fls = withTempDir dir $ \tmpD -> do
      let sourceCells = getZipSource $ (,,) <$>
              ZipSource (iterateC succ 0) <*>
              ZipSource seqDepthC <*>
              ZipSource ( sourceFile (tmpD <> "/embed") .| linesUnboundedAsciiC .|
                mapC (map readDouble . B.split '\t') )
          input = T.pack $ intercalate "," $ map (^.location) mats
      shelly $ run_ "sc_utils" $ [ "clust", input, T.pack tmpD <> "/clust",
          "--embed", T.pack tmpD <> "/embed" ] ++ toParams opt
      cells <- runResourceT $ runConduit $ sourceCells .| mapC f .| sinkVector
      clusters <- readClusters $ tmpD <> "/clust"
      return $ zipWith (\i -> CellCluster $ B.pack $ "C" ++ show i) [1::Int ..] $
          map (map (cells V.!)) clusters
  where
    coverage = head $ fst $ unzip fls
    mats = snd $ unzip fls
    readClusters fl = map (map readInt . B.split ',') . B.lines <$>
        B.readFile fl
    seqDepthC = sourceFile (coverage^.location) .| linesUnboundedAsciiC .|
        mapC ((\[a,b] -> (a,b)) . B.split '\t')
    f (i, (bc, dep), [d1,d2,d3,d4,d5]) = Cell i (d1,d2) (d3,d4,d5) bc $ readInt dep
    f _ = error "formatting error"
{-# INLINE clustering' #-}

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
                S.toList $ S.fromList $ M.elems clusters
        fileHandles <- mapM (\x -> openFile x WriteMode) outputs
        let f x = let h = M.lookupDefault undefined bc fileHandles
                      bc = M.lookupDefault (error $ show nm) nm clusters
                      nm = fromJust ((x::BED)^.name) 
                  in liftIO $ B.hPutStrLn h $ toLine x 
        runResourceT $ runConduit $ streamBedGzip (bed^.location) .| mapM_C f
        return $ M.toList $ fmap (\x -> location .~ x $ emptyFile) outputs )
  where
    getClusterBarcodes :: B.ByteString    -- ^ experiment id
                    -> SCATACSeq S (File '[] 'Other)   -- ^ Cluster files
                    -> IO (M.HashMap B.ByteString B.ByteString)  -- ^ Barcode to cluster map
    getClusterBarcodes exp_id clFl = do
        clusters <- decodeFile $ clFl^.replicates._2.files.location :: IO [CellCluster]
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
        , plotDupRate res ]
      where
        res = flip map cls $ \x ->
            (T.pack $ B.unpack $ _cluster_name x, map f $ _cluster_member x)
        f x = M.lookupDefault undefined (_cell_barcode x) statMap
        statMap = M.fromList $ map (\x -> (_barcode x, x)) stats

visualizeCluster :: [CellCluster] -> [EChart]
visualizeCluster cs = [scatter' dat2D <> toolbox, scatter' dat2D' <> toolbox] 
  where
    dat2D = flip map cs $ \(CellCluster nm cells) ->
        (B.unpack nm, map _cell_2d cells)
    dat2D' = map (first head . unzip) $ groupBy ((==) `on` fst) $ sortBy (comparing fst) $ concatMap
        (map (\x -> (getName $ _cell_barcode x, _cell_2d x)) . _cluster_member) cs
    getName x = let prefix = fst $ B.breakEnd (=='+') x
                in if B.null prefix then "" else B.unpack $ B.init prefix

-- | Compute the normalized tissue composition for each cluster.
tissueComposition :: DF.DataFrame Int -> EChart
tissueComposition = stackBar . DF.map round' . DF.mapCols normalize .
    DF.transpose . DF.mapCols normalize . DF.map fromIntegral
  where
    round' x = fromIntegral (round $ x * 1000 :: Int) / 1000
    normalize xs | V.all (==0) xs = xs
                 | otherwise = V.map (/V.sum xs) xs

-- | Compute the cluster composition for each tissue.
clusterComposition :: DF.DataFrame Int -> EChart
clusterComposition = stackBar . DF.map round' . DF.mapCols normalize .
    DF.map fromIntegral
  where
    round' x = fromIntegral (round $ x * 1000 :: Int) / 1000
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

-- | Random sample 30,000 cells.
sampleCells :: [CellCluster] -> IO [CellCluster]
sampleCells clusters
    | ratio >= 1 = return clusters
    | otherwise = do
        gen <- create
        forM clusters $ \c -> do
            s <- sampling gen ratio $ V.fromList $ _cluster_member c
            return $ c {_cluster_member = V.toList s}
  where
    n = foldl1' (+) $ map (length . _cluster_member) clusters
    ratio = 1 / (fromIntegral n / 30000) :: Double
    sampling gen frac v = V.take n' <$> uniformShuffle v gen
      where
        n' = max 200 $ truncate $ frac * fromIntegral (V.length v)