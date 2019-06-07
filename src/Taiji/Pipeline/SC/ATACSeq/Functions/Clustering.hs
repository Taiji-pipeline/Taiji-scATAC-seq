{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE DataKinds #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Clustering
    ( module Taiji.Pipeline.SC.ATACSeq.Functions.Clustering.LSA
    , module Taiji.Pipeline.SC.ATACSeq.Functions.Clustering.LDA
    , plotClusters
    , plotClusters'
    , clust
    ) where

import qualified Data.ByteString.Char8 as B
import qualified Data.Text as T
import qualified Data.HashMap.Strict as M
import qualified Data.Vector as V
import System.Random.MWC.Distributions
import System.Random.MWC
import Bio.Utils.Misc (readDouble, readInt)
import Shelly (shelly, run_)
   
import Taiji.Pipeline.SC.ATACSeq.Functions.Clustering.LSA
import Taiji.Pipeline.SC.ATACSeq.Functions.Clustering.LDA
import Taiji.Prelude
import Taiji.Pipeline.SC.ATACSeq.Types
import Taiji.Pipeline.SC.ATACSeq.Functions.Utils

clust :: Maybe FilePath      -- ^ temp dir
      -> (File '[] 'Tsv, File '[Gzip] 'Tsv)   -- ^ lsa input matrix
      -> IO [CellCluster]
clust dir (coverage, mat) = withTempDir dir $ \tmpD -> do
      runResourceT $ runConduit $ seqDepthC .| mapC snd .|
          unlinesAsciiC .| sinkFile (tmpD <> "/coverage")
      shelly $ run_ "sc_utils" [ "clust",
          "--discard",
          --"--coverage", T.pack tmpD <> "/coverage",
          "--embed", T.pack tmpD <> "/embed",
          T.pack $ mat^.location, T.pack tmpD <> "/clust" ]

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

plotClusters :: FilePath
             -> SCATACSeq S [CellCluster]
             -> IO ()
plotClusters dir input = do
    let output = printf "%s/%s_rep%d_cluster.html" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    visualizeCluster output Nothing $ input^.replicates._2.files

plotClusters' :: SCATACSeqConfig config
              => [CellCluster]
              -> ReaderT config IO ()
plotClusters' clusters = do
    dir <- asks ((<> "/Cluster") . _scatacseq_output_dir) >>= getPath
    let output = dir <> "/cluster.html"
        bcs = Just $ V.fromList $ concatMap
            (map (B.unpack . fst . B.break (=='_') . _cell_barcode) . _cluster_member) clusters
    liftIO $ sampleCells clusters >>= visualizeCluster output bcs

sampleCells :: [CellCluster]
            -> IO [CellCluster]
sampleCells clusters = do
    gen <- create
    forM clusters $ \c -> do
        s <- sampling gen 0.2 $ V.fromList $ _cluster_member c
        return $ c {_cluster_member = V.toList s}

sampling :: GenIO
         -> Double  -- ^ fraction
         -> V.Vector a -> IO (V.Vector a)
sampling gen frac v = V.take n <$> uniformShuffle v gen
  where
    n = truncate $ frac * fromIntegral (V.length v)