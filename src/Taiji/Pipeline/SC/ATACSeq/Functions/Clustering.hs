{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE DeriveLift #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Clustering
    ( plotClusters
    , clust
    , lsaBuilder
    , dmClust
    , extractTags
    , extractSubMatrix
    , doClustering
    , getBedCluster
    , getClusterBarcodes
    , extractBedByBarcode 
    , Embedding(..)
    , Normalization(..)
    , ClustOpt(..)
    , defClustOpt
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
import Bio.Utils.Misc (readDouble, readInt)
import Shelly (shelly, run_, escaping)
import Control.Workflow
import Data.Conduit.Zlib (gzip)
   
import Taiji.Pipeline.SC.ATACSeq.Functions.Feature
import Taiji.Pipeline.SC.ATACSeq.Functions.DR.LSA
import Taiji.Pipeline.SC.ATACSeq.Functions.DR.DiffusionMap
import Taiji.Pipeline.SC.ATACSeq.Functions.Clustering.Utils
import Taiji.Prelude
import Taiji.Pipeline.SC.ATACSeq.Types
import Taiji.Pipeline.SC.ATACSeq.Functions.Utils

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
defClustOpt = ClustOpt UnitBall UMAP Nothing 20 Nothing

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

 
doClustering :: SCATACSeqConfig config
             => FilePath
             -> ClustOpt
             -> SCATACSeq S (File '[] 'Tsv, File '[Gzip] 'Tsv)
             -> ReaderT config IO (SCATACSeq S (File '[] 'Other))
doClustering prefix opt input = do
    tmp <- asks _scatacseq_temp_dir
    dir <- asks ((<> asDir ("/" ++ prefix)) . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_clusters.bin" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traversed.files %%~ liftIO . ( \fl -> do
        clust opt tmp fl >>= encodeFile output
        return $ location .~ output $ emptyFile )

clust :: ClustOpt
      -> Maybe FilePath      -- ^ temp dir
      -> (File '[] 'Tsv, File '[Gzip] 'Tsv)   -- ^ lsa input matrix
      -> IO [CellCluster]
clust opt dir (coverage, mat) = withTempDir dir $ \tmpD -> do
      let sourceCells = getZipSource $ (,,) <$>
              ZipSource (iterateC succ 0) <*>
              ZipSource seqDepthC <*>
              ZipSource ( sourceFile (tmpD <> "/embed") .| linesUnboundedAsciiC .|
                mapC (map readDouble . B.split '\t') )
      shelly $ run_ "sc_utils" $ [ "clust",
          T.pack $ mat^.location, T.pack tmpD <> "/clust",
          "--embed", T.pack tmpD <> "/embed" ] ++ toParams opt
      cells <- runResourceT $ runConduit $ sourceCells .| mapC f .| sinkVector
      clusters <- readClusters $ tmpD <> "/clust"
      return $ zipWith (\i -> CellCluster $ B.pack $ "C" ++ show i) [1::Int ..] $
          map (map (cells V.!)) clusters
  where
    readClusters fl = map (map readInt . B.split ',') . B.lines <$>
        B.readFile fl
    seqDepthC = sourceFile (coverage^.location) .| linesUnboundedAsciiC .|
        mapC ((\[a,b] -> (a,b)) . B.split '\t')
    f (i, (bc, dep), [d1,d2,d3,d4,d5]) = Cell i (d1,d2) (d3,d4,d5) bc $ readInt dep
    f _ = error "formatting error"
{-# INLINE clust #-}

--------------------------------------------------------------------------------
-- LSA
--------------------------------------------------------------------------------
lsaBuilder :: Builder ()
lsaBuilder = do
    -- Clustering 1st round (By Window)
    namespace "Window" $ lsaClust "/Cluster_by_window/LSA/" defClustOpt
    path ["Merge_Window_Matrix", "Window_LSA_Reduce"]

    -- Extract tags for each cluster
    namespace "Window_LSA" $ extractTags "/temp/Bed/Cluster/"
    ["Get_Bed", "Window_LSA_Cluster"] ~> "Window_LSA_Extract_Tags_Prep"

    -- Call peaks 1st round
    genPeakMat "/temp/Peak/" (Just "LSA_1st") 
        "Window_LSA_Merge_Tags" "Get_Windows"

    -- Clustering 2nd round
    namespace "Peak" $ lsaClust "/Cluster_by_peak/LSA/" defClustOpt
    path ["LSA_1st_Merge_Peak_Matrix", "Peak_LSA_Reduce"]

    -- Subclustering
    node "Extract_Sub_Matrix" [| \(x,y) -> 
        let [input] = zipExp x y
        in extractSubMatrix "/temp/" input |] $ return ()
    ["LSA_1st_Merge_Peak_Matrix", "Peak_LSA_Cluster"] ~> "Extract_Sub_Matrix"
    namespace "SubCluster" $ lsaClust "/Cluster_by_peak/LSA/SubCluster/" $ defClustOpt{_resolution = Just 0.0007}
    path ["Extract_Sub_Matrix", "SubCluster_LSA_Reduce"]


-- | Perform LSA analysis.
lsaClust :: FilePath   -- ^ Directory to save the results
         -> ClustOpt
         -> Builder ()
lsaClust prefix opt = do
    nodePar "LSA_Reduce" [| performLSA prefix |] $ return ()
    nodePar "LSA_Cluster" [| \x -> do
        res <- asks _scatacseq_cluster_resolution 
        case _resolution opt of
            Nothing -> doClustering prefix opt{_resolution=res} x
            _ -> doClustering prefix opt x
        |] $ return ()
    nodePar "LSA_Viz" [| \x -> do
        dir <- asks ((<> asDir ("/" ++ prefix)) . _scatacseq_output_dir) >>= getPath
        liftIO $ plotClusters dir x
        |] $ return ()
    path ["LSA_Reduce", "LSA_Cluster", "LSA_Viz"]


-- | Perform LSA analysis.
dmClust :: FilePath   -- ^ Directory to save the results
        -> ClustOpt
        -> Builder ()
dmClust prefix opt = do
    nodePar "DM_Reduce" [| performDM prefix |] $ return ()
    nodePar "DM_Cluster" [| doClustering prefix opt |] $ return ()
    nodePar "DM_Viz" [| \x -> do
        dir <- asks ((<> asDir ("/" ++ prefix)) . _scatacseq_output_dir) >>= getPath
        liftIO $ plotClusters dir x
        |] $ return ()
    path ["DM_Reduce", "DM_Cluster", "DM_Viz"]

-- | Extract tags for clusters.
extractTags :: FilePath   -- ^ Directory to save the results
            -> Builder ()
extractTags prefix = do
    node "Extract_Tags_Prep"  [| liftIO . getClusterBarcodes |] $ return ()
    nodePar "Extract_Tags" 'getBedCluster $ return ()
    node "Merge_Tags_Prep" [| \input -> return $
        map (first head . unzip) $ groupBy ((==) `on` fst) $
        sortBy (comparing fst) $ concat $ input^..folded.replicates._2.files
        |] $ return ()
    nodePar "Merge_Tags" [| mergeBedCluster prefix |] $ return ()
    path ["Extract_Tags_Prep", "Extract_Tags", "Merge_Tags_Prep", "Merge_Tags"]

-- | Get barcodes
getClusterBarcodes :: ([SCATACSeq S file], [SCATACSeq S (File '[] 'Other)])
                   -> IO [(SCATACSeq S (file, [(B.ByteString, [B.ByteString])]))]
getClusterBarcodes (inputs, clusters) = do
    clusters' <- case clusters of
        [x] -> zip (repeat "") <$> decodeFile (x^.replicates._2.files.location)
        xs -> fmap concat $ forM xs $ \x ->
            zip (repeat $ B.pack $ T.unpack $ x^.eid) <$>
            decodeFile (x^.replicates._2.files.location)
    return $ flip map inputs $ \input ->
        let res = flip map clusters' $ \(prefix, CellCluster{..}) ->
                (prefix <> _cluster_name, mapMaybe (getBarcode (input^.eid)) _cluster_member)
        in input & replicates.traverse.files %~ (\f -> (f, res))
  where
    getBarcode e x | B.unpack (B.init i) == T.unpack e = Just bc
                   | otherwise = Nothing
      where
        (i, bc) = B.breakEnd (=='+') $ _cell_barcode x

-- | Extract BEDs for each cluster.
getBedCluster :: SCATACSeqConfig config
              => SCATACSeq S ( File '[NameSorted, Gzip] 'Bed
                             , [(B.ByteString, [B.ByteString])] )  -- ^ clusters
              -> ReaderT config IO
                    (SCATACSeq S [(B.ByteString, File '[] 'Bed)])
getBedCluster input = do
    let idRep = asDir $ "/temp/" <> T.unpack (input^.eid) <>
            "_rep" <> show (input^.replicates._1)
    dir <- asks _scatacseq_output_dir >>= getPath . (<> idRep)
    input & replicates.traverse.files %%~ liftIO . ( \(bed, _) -> do
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
    clusters = M.fromListWith (error "same barcode") $
        concatMap (\(x, ys) -> zip ys $ repeat x) $ input^.replicates._2.files._2

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

plotClusters :: FilePath
             -> SCATACSeq S (File '[] 'Other)
             -> IO ()
plotClusters dir input = do
    inputData <- decodeFile $ input^.replicates._2.files.location
    let output = printf "%s/%s_rep%d_cluster.html" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
        stat = printf "%s/%s_rep%d_cluster_stat.html" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    when (B.elem '+' $ _cell_barcode $ head $ _cluster_member $ head inputData) $ clusterStat stat inputData
    clusters <- sampleCells inputData
    visualizeCluster output clusters

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
        n' = max 100 $ truncate $ frac * fromIntegral (V.length v)

-- | Extract submatrix
extractSubMatrix :: SCATACSeqConfig config
                 => FilePath   -- ^ dir
                 -> SCATACSeq S ((a, File tags 'Other), File '[] 'Other)
                 -> ReaderT config IO [SCATACSeq S (a, File tags 'Other)]
extractSubMatrix prefix input = do
    dir <- asks _scatacseq_output_dir >>= getPath . (<> (asDir prefix))
    let ((idxFl, matFl), clFl) = input^.replicates._2.files
    liftIO $ do
        clusters <- decodeFile $ clFl^.location
        mat <- mkSpMatrix id $ matFl^.location
        let mkSink CellCluster{..} = filterC ((`S.member` ids) . fst) .|
                (yield header >> mapC (encodeRowWith id)) .| unlinesAsciiC .|
                gzip .| (sinkFile output >> return res)
              where
                res = input & eid .~ input^.eid <> "+" <> T.pack (B.unpack _cluster_name)
                            & replicates._2.files .~ (idxFl, location .~ output $ matFl)
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