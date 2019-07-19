{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE DataKinds #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Clustering
    ( module Taiji.Pipeline.SC.ATACSeq.Functions.Clustering.LSA
    , module Taiji.Pipeline.SC.ATACSeq.Functions.Clustering.LDA
    , module Taiji.Pipeline.SC.ATACSeq.Functions.Clustering.SnapTools
    , module Taiji.Pipeline.SC.ATACSeq.Functions.Clustering.DiffusionMap
    , plotClusters
    , clust
    , lsaClust
    , ldaClust
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
    ) where

import qualified Data.ByteString.Char8 as B
import Data.Binary (encodeFile, decodeFile)
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
   
import Taiji.Pipeline.SC.ATACSeq.Functions.Clustering.LSA
import Taiji.Pipeline.SC.ATACSeq.Functions.Clustering.LDA
import Taiji.Pipeline.SC.ATACSeq.Functions.Clustering.SnapTools
import Taiji.Pipeline.SC.ATACSeq.Functions.Clustering.DiffusionMap
import Taiji.Pipeline.SC.ATACSeq.Functions.Clustering.Utils
import Taiji.Prelude
import Taiji.Pipeline.SC.ATACSeq.Types
import Taiji.Pipeline.SC.ATACSeq.Functions.Utils

-- | Perform LSA analysis.
lsaClust :: FilePath   -- ^ Directory to save the results
         -> Builder ()
lsaClust prefix = do
    nodePar "LSA_Reduce" [| performLSA prefix |] $ return ()
    nodePar "LSA_Cluster" [| doClustering prefix $ ClustOpt UnitBall UMAP |] $ return ()
    nodePar "LSA_Viz" [| \x -> do
        dir <- asks ((<> asDir ("/" ++ prefix)) . _scatacseq_output_dir) >>= getPath
        liftIO $ plotClusters dir x
        |] $ return ()
    path ["LSA_Reduce", "LSA_Cluster", "LSA_Viz"]

-- | Perform LDA analysis.
ldaClust :: FilePath   -- ^ Directory to save the results
         -> Builder ()
ldaClust prefix = do
    nodePar "LDA_Reduce" [| performLDA prefix |] $ return ()
    nodePar "LDA_Cluster" [| doClustering prefix $ ClustOpt None UMAP |] $ return ()
    nodePar "LDA_Viz" [| \x -> do
        dir <- asks ((<> asDir ("/" ++ prefix)) . _scatacseq_output_dir) >>= getPath
        liftIO $ plotClusters dir x
        |] $ return ()
    path ["LDA_Reduce", "LDA_Cluster", "LDA_Viz"]

-- | Perform LSA analysis.
dmClust :: FilePath   -- ^ Directory to save the results
        -> Builder ()
dmClust prefix = do
    nodePar "DM_Reduce" [| performDM prefix |] $ return ()
    nodePar "DM_Cluster" [| doClustering prefix $ ClustOpt None UMAP
        |] $ return ()
    nodePar "DM_Viz" [| \x -> do
        dir <- asks ((<> asDir ("/" ++ prefix)) . _scatacseq_output_dir) >>= getPath
        liftIO $ plotClusters dir x
        |] $ return ()
    path ["DM_Reduce", "DM_Cluster", "DM_Viz"]
  where

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

data Embedding = UMAP
               | TSNE

data Normalization = Drop1st
                   | UnitBall
                   | Regression
                   | None

data ClustOpt = ClustOpt
    { _normalization :: Normalization
    , _embedding_method :: Embedding
    }

clust :: ClustOpt
      -> Maybe FilePath      -- ^ temp dir
      -> (File '[] 'Tsv, File '[Gzip] 'Tsv)   -- ^ lsa input matrix
      -> IO [CellCluster]
clust ClustOpt{..} dir (coverage, mat) = withTempDir dir $ \tmpD -> do
      runResourceT $ runConduit $ seqDepthC .| mapC snd .|
          unlinesAsciiC .| sinkFile (tmpD <> "/coverage")
      let embed = case _embedding_method of
              UMAP -> ["--embed-method", "umap"]
              TSNE -> ["--embed-method", "tsne"]
          normalize = case _normalization of
              Drop1st -> ["--discard"]
              Regression -> ["--coverage", T.pack tmpD <> "/coverage"]
              UnitBall -> ["--scale"]
              None -> []
          sourceCells = getZipSource $ (,,) <$>
              ZipSource (iterateC succ 0) <*>
              ZipSource seqDepthC <*>
              ZipSource ( sourceFile (tmpD <> "/embed") .| linesUnboundedAsciiC .|
                mapC (map readDouble . B.split '\t') )
      shelly $ run_ "sc_utils" $ [ "clust", T.pack $ mat^.location
          , T.pack tmpD <> "/clust", "--embed", T.pack tmpD <> "/embed" ] ++
          embed ++ normalize
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

extractBedByBarcode :: Elem 'Gzip tags ~ 'True
                    => [FilePath]   -- ^ Outputs
                    -> [[B.ByteString]]  -- ^ Barcodes
                    -> File tags 'Bed
                    -> IO [File '[] 'Bed]
extractBedByBarcode outputs bcs input = do
    fileHandles <- V.fromList <$> mapM (\x -> openFile x WriteMode) outputs
    let bcIdx = M.fromList $ concat $ zipWith (\i x -> zip x $ repeat i) [0..] bcs
        f x = let idx = M.lookupDefault (error $ show nm) nm bcIdx
                  nm = fromJust ((x::BED)^.name) 
              in liftIO $ B.hPutStrLn (fileHandles V.! idx) $ toLine x 
    runResourceT $ runConduit $ streamBedGzip (input^.location) .| mapM_C f
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

sampleCells :: [CellCluster]
            -> IO [CellCluster]
sampleCells clusters
    | ratio >= 1 = return clusters
    | otherwise = do
        gen <- create
        forM clusters $ \c -> do
            s <- sampling gen ratio $ V.fromList $ _cluster_member c
            return $ c {_cluster_member = V.toList s}
  where
    n = foldl1' (+) $ map (length . _cluster_member) clusters
    ratio = 1 / (fromIntegral n / 30000)

sampling :: GenIO
         -> Double  -- ^ fraction
         -> V.Vector a -> IO (V.Vector a)
sampling gen frac v = V.take n <$> uniformShuffle v gen
  where
    n = truncate $ frac * fromIntegral (V.length v)

extractSubMatrix :: SCATACSeqConfig config
                 => ( [SCATACSeq S (File tags 'Other)]
                    , [SCATACSeq S (File '[] 'Other)] )
                 -> ReaderT config IO [SCATACSeq S (File tags 'Other)]
extractSubMatrix ([input], [clFl]) = do
    dir <- asks _scatacseq_output_dir >>= getPath . (<> "/Feature/Peak/Cluster/")
    liftIO $ do
        clusters <- decodeFile $ clFl^.replicates._2.files.location
        mat <- mkSpMatrix id $ input^.replicates._2.files.location
        let mkSink CellCluster{..} = filterC ((`S.member` ids) . fst) .|
                (yield header >> mapC (encodeRowWith id)) .| unlinesAsciiC .|
                gzip .| (sinkFile output >> return res)
              where
                res = input & eid .~ T.pack (B.unpack _cluster_name)
                            & replicates._2.files.location .~ output
                output = dir <> B.unpack _cluster_name <> "_subMat.txt.gz"
                header = B.pack $ printf "Sparse matrix: %d x %d" nCell
                    (_num_col mat)
                ids = S.fromList $ map _cell_barcode _cluster_member
                nCell = S.size ids
        runResourceT $ runConduit $ streamRows mat .|
            sequenceSinks (map mkSink clusters)
extractSubMatrix _ = return []