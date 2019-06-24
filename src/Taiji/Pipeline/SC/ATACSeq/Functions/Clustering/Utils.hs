{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE QuasiQuotes       #-}
{-# LANGUAGE ViewPatterns #-}
{-# LANGUAGE GADTs #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Clustering.Utils
    ( clusterStat
    ) where

import Data.Int (Int32)
import qualified Data.HashMap.Strict as M
import qualified Data.ByteString.Char8 as B
import qualified Data.Vector as V
import Data.List
import Data.List.Ordered (nubSort)
import Data.Ord
import Data.Function (on)
import qualified Data.Text as T
{-
import qualified Language.R                        as R
import           Language.R.QQ
import qualified Data.Vector.SEXP as V
import Language.R.HExp
-}

import Taiji.Pipeline.SC.ATACSeq.Functions.Utils
import qualified Taiji.Utils.DataFrame as DF
import Taiji.Utils.Plot.ECharts
import Taiji.Utils.Plot
import Taiji.Prelude

{-
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
    -}

clusterStat :: FilePath -> [CellCluster] -> IO ()
clusterStat output clusters = savePlots output []
    [ stackBar $ DF.mapCols normalize df
    , stackBar $ DF.mapCols normalize $ DF.transpose df
    , heatmap $ DF.orderDataFrame id $ DF.spearman df
    , heatmap $ DF.orderDataFrame id $ DF.spearman $ DF.transpose df ]
  where
    df = DF.mkDataFrame rownames colnames $
        map (\x -> map (\i -> fromIntegral $ M.lookupDefault 0 i x) colnames) rows
    (rownames, rows) = unzip $ map f clusters
    colnames = nubSort $ concatMap M.keys rows
    f CellCluster{..} = ( T.pack $ B.unpack _cluster_name,
        M.fromListWith (+) $ map (\x -> (tissueName x, 1)) _cluster_member )
    tissueName Cell{..} = T.pack $ B.unpack $ B.init $ fst $
        B.breakEnd (=='+') _cell_barcode
    normalize xs = V.map (\x -> round' $ x / s) xs
      where 
        round' x = fromIntegral (round $ x * 1000) / 1000
        s = V.foldl1' (+) xs