{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE QuasiQuotes       #-}
{-# LANGUAGE ViewPatterns #-}
{-# LANGUAGE GADTs #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Clustering.Utils
    ( seurat
    , clusterStat
    ) where

import Data.Int (Int32)
import qualified Data.HashMap.Strict as M
import qualified Data.ByteString.Char8 as B
import Data.List
import Data.Ord
import Data.Function (on)
import qualified Language.R                        as R
import           Language.R.QQ
import qualified Data.Vector.SEXP as V
import Language.R.HExp

import Taiji.Pipeline.SC.ATACSeq.Functions.Utils
import Taiji.Utils.Plot.ECharts
import Taiji.Utils.Plot
import Taiji.Prelude

-- | Perform clustering using the Seurat library (SNN+graph_modularity).
seurat :: [String]   -- ^ Sample barcodes/names.
       -> [[Double]]   -- ^ Columns are features, rows are samples.
       -> IO [[String]]
seurat names xs = R.runRegion $ do
    membership <- [r| library("Seurat")
        input <- matrix(xs__hs, nrow=ncol_hs, byrow=T)
        colnames(input) <- names_hs
        rownames(input) <- features_hs
        seurat_obj <- Seurat::CreateSeuratObject(input)
        seurat_obj[['pca']] <- Seurat::CreateDimReducObject(embeddings=t(input), key='PC_', assay='RNA')
        seurat_obj <- Seurat::FindNeighbors(seurat_obj, reduction="pca", k.param=3,
            nn.eps=0, dims=1:ncol_hs, do.plot=T, graph.name="grx")
        print(seurat_obj$grx)
        seurat_obj <- Seurat::FindClusters(seurat_obj, n.start=20,
            resolution=0.3, dims=1:ncol_hs, reduction="pca")
        seurat_obj$seurat_clusters
    |]
    let members = R.unSomeSEXP membership $
            \(hexp -> Int vec) -> V.toList vec
    return $ map (fst . unzip) $ groupBy ((==) `on` snd) $
        sortBy (comparing snd) $ zip names (members :: [Int32])
  where
    xs_ = concat xs
    features = map show $ [1 .. length (head xs)]
    ncol = fromIntegral $ length $ head xs :: Double

clusterStat :: FilePath -> [CellCluster] -> IO ()
clusterStat output clusters = savePlots output [] [chart]
  where
    chart = stackBar "" (map B.unpack allTissues) $
        map (\(c, x) -> (B.unpack c, map (\t -> fromIntegral $ M.lookupDefault 0 t x) allTissues)) byCluster
    allTissues = map fst byTissue
    allClusters = map fst byCluster
    byCluster :: [(B.ByteString, M.HashMap B.ByteString Int)]
    byCluster = M.toList $ foldl' f M.empty clusters
      where
        f m CellCluster{..} = foldl' g m $ map (\x -> (_cluster_name, tissueName x)) _cluster_member
        g m (k, v) = M.alter fun k m
          where
            fun Nothing = Just $ M.singleton v 1
            fun (Just x) = Just $ M.insertWith (+) v 1 x
    byTissue :: [(B.ByteString, M.HashMap B.ByteString Int)]
    byTissue = M.toList $ foldl' f M.empty clusters
      where
        f m CellCluster{..} = foldl' g m $ map (\x -> (tissueName x, _cluster_name)) _cluster_member
        g m (k, v) = M.alter fun k m
          where
            fun Nothing = Just $ M.singleton v 1
            fun (Just x) = Just $ M.insertWith (+) v 1 x
    tissueName Cell{..} = B.init $ fst $ B.breakEnd (=='_') _cell_barcode