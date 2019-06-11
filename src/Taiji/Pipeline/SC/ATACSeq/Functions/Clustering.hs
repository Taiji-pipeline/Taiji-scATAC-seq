{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE DataKinds #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Clustering
    ( module Taiji.Pipeline.SC.ATACSeq.Functions.Clustering.LSA
    , module Taiji.Pipeline.SC.ATACSeq.Functions.Clustering.LDA
    , plotClusters
    , clust
    , lsaClust
    , extractTags
    ) where

import qualified Data.ByteString.Char8 as B
import qualified Data.Text as T
import qualified Data.HashSet as S
import Bio.Data.Bed
import Control.Arrow (first)
import qualified Data.Vector as V
import System.Random.MWC.Distributions
import System.Random.MWC
import Bio.Utils.Misc (readDouble, readInt)
import Shelly (shelly, run_, escaping)
import Control.Workflow
   
import Taiji.Pipeline.SC.ATACSeq.Functions.Clustering.LSA
import Taiji.Pipeline.SC.ATACSeq.Functions.Clustering.LDA
import Taiji.Pipeline.SC.ATACSeq.Functions.Clustering.Utils
import Taiji.Prelude
import Taiji.Pipeline.SC.ATACSeq.Types
import Taiji.Pipeline.SC.ATACSeq.Functions.Utils

lsaClust :: FilePath -> Builder ()
lsaClust prefix = do
    nodePar "LSA" [| performLSA prefix |] $ return ()
    nodePar "Cluster_LSA" [| \input -> do
        tmp <- asks _scatacseq_temp_dir
        input & replicates.traversed.files %%~ liftIO . clust True tmp
        |] $ return ()
    nodePar "Visualize_LSA_Cluster" [| \x -> do
        dir <- asks ((<> asDir ("/" ++ prefix)) . _scatacseq_output_dir) >>= getPath
        liftIO $ plotClusters dir x
        |] $ return ()
    path ["LSA", "Cluster_LSA", "Visualize_LSA_Cluster"]

extractTags :: Builder ()
extractTags = do
    node "Extract_Tags_Prep"  [| return . getClusterBarcodes |] $ return ()
    nodePar "Extract_Tags" 'getBedCluster $ return ()
    node "Merge_Tags_Cluster" 'mergeBedCluster $ return ()
    path ["Extract_Tags_Prep", "Extract_Tags", "Merge_Tags_Cluster"]

clust :: Bool   -- ^ Whether to discard the first dimension
      -> Maybe FilePath      -- ^ temp dir
      -> (File '[] 'Tsv, File '[Gzip] 'Tsv)   -- ^ lsa input matrix
      -> IO [CellCluster]
clust discard dir (coverage, mat) = withTempDir dir $ \tmpD -> do
      runResourceT $ runConduit $ seqDepthC .| mapC snd .|
          unlinesAsciiC .| sinkFile (tmpD <> "/coverage")
      shelly $ run_ "sc_utils" $ [ "clust",
          --"--coverage", T.pack tmpD <> "/coverage",
          "--embed", T.pack tmpD <> "/embed",
          T.pack $ mat^.location, T.pack tmpD <> "/clust" ] ++
            if discard then ["--discard"] else []

      let sourceCells = getZipSource $ (,,) <$>
              ZipSource (iterateC succ 0) <*>
              ZipSource seqDepthC <*>
              ZipSource ( sourceFile (tmpD <> "/embed") .| linesUnboundedAsciiC .|
                mapC (map readDouble . B.split '\t') )
      cells <- runResourceT $ runConduit $ sourceCells .| mapC f .| sinkVector
      clusters <- readClusters $ tmpD <> "/clust"
      return $ zipWith (\i -> CellCluster $ B.pack $ "C" ++ show i) [1::Int ..] $
          map (map (cells V.!)) clusters
  where
    readClusters fl = map (map readInt . B.split ',') . B.lines <$>
        B.readFile fl
    seqDepthC = sourceFile (coverage^.location) .| linesUnboundedAsciiC .|
        mapC ((\[a,b] -> (a,b)) . B.split '\t')
    f (i, (bc, dep), [x,y,z]) = Cell i x y z bc $ readInt dep
    f _ = error "formatting error"
{-# INLINE clust #-}

getClusterBarcodes :: ([SCATACSeq S file], [SCATACSeq S [CellCluster]])
                   -> [(SCATACSeq S (file, [(B.ByteString, [B.ByteString])]))]
getClusterBarcodes (inputs, [clusters]) = flip map inputs $ \input ->
    let res = flip map (clusters^.replicates._2.files) $ \CellCluster{..} ->
            (_cluster_name, mapMaybe (getBarcode (input^.eid)) _cluster_member)
    in input & replicates.traverse.files %~ (\f -> (f, res))
  where
    getBarcode e x | B.unpack (B.init i) == T.unpack e = Just bc
                   | otherwise = Nothing
      where
        (i, bc) = B.breakEnd (=='_') $ _cell_barcode x
getClusterBarcodes _ = []

-- | Extract BEDs for each cluster.
getBedCluster :: SCATACSeqConfig config
              => SCATACSeq S ( File '[NameSorted, Gzip] 'Bed
                             , [(B.ByteString, [B.ByteString])] )  -- ^ clusters
              -> ReaderT config IO
                    (SCATACSeq S [(B.ByteString, File '[Gzip] 'Bed)])
getBedCluster input = do
    let idRep = asDir $ "/temp/" <> T.unpack (input^.eid) <>
            "_rep" <> show (input^.replicates._1)
    dir <- asks _scatacseq_output_dir >>= getPath . (<> idRep)
    input & replicates.traverse.files %%~ liftIO . ( \(bed, cs) -> do
        let sinks = sequenceConduits $ flip map cs $ \(cName, bc) -> do
                let output = dir ++ "/" ++ B.unpack cName ++ ".bed.gz"
                    cells = S.fromList bc
                    fl = location .~ output $ emptyFile
                filterC (\x -> fromJust ((x::BED)^.name) `S.member` cells) .|
                    sinkFileBedGzip output
                return (cName, fl)
        runResourceT $ runConduit $
            streamBedGzip (bed^.location) .| sinks )

-- | Extract BEDs for each cluster.
mergeBedCluster :: SCATACSeqConfig config
                => [SCATACSeq S [(B.ByteString, File '[Gzip] 'Bed)]]
                -> ReaderT config IO [(B.ByteString, File '[Gzip] 'Bed)]
mergeBedCluster inputs = do
    dir <- asks _scatacseq_output_dir >>= getPath . (<> "/Bed/Cluster/")
    liftIO $ forM clusters $ \(cName, fls) -> do
        let output = dir <> B.unpack cName <> ".bed.gz"
        shelly $ escaping False $ do
            run_ "cat" $ map (\x -> T.pack $ x^.location) fls ++ [">", T.pack output]
            run_ "rm" $ map (\x -> T.pack $ x^.location) fls
        return (cName, location .~ output $ emptyFile)
  where
    clusters = map (first head . unzip) $ groupBy ((==) `on` fst) $
        sortBy (comparing fst) $ concat $ inputs^..folded.replicates._2.files

plotClusters :: FilePath
             -> SCATACSeq S [CellCluster]
             -> IO ()
plotClusters dir input = do
    let output2d = printf "%s/%s_rep%d_cluster_2d.html" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
        output3d = printf "%s/%s_rep%d_cluster_3d.html" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
        stat = printf "%s/%s_rep%d_cluster_stat.html" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    clusterStat stat inputData
    clusters <- sampleCells inputData
    visualizeCluster2D output2d Nothing clusters
    visualizeCluster3D output3d Nothing clusters
  where
    inputData = input^.replicates._2.files

{-
plotClusters' :: SCATACSeqConfig config
              => [CellCluster]
              -> ReaderT config IO ()
plotClusters' clusters = do
    dir <- asks ((<> "/Cluster") . _scatacseq_output_dir) >>= getPath
    let output = dir <> "/cluster.html"
        bcs = Just $ V.fromList $ concatMap
            (map (B.unpack . fst . B.break (=='_') . _cell_barcode) . _cluster_member) clusters
    liftIO $ sampleCells clusters >>= visualizeCluster output bcs
    -}

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