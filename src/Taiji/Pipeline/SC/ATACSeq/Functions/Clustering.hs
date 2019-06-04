{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE DataKinds #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Clustering
    ( module Taiji.Pipeline.SC.ATACSeq.Functions.Clustering.LSA
    , plotClusters
    ) where

import qualified Data.ByteString.Char8 as B
import qualified Data.Text as T
import qualified Data.HashMap.Strict as M
   
import Taiji.Pipeline.SC.ATACSeq.Functions.Clustering.LSA
import Taiji.Prelude
import Taiji.Pipeline.SC.ATACSeq.Types
import Taiji.Pipeline.SC.ATACSeq.Functions.Utils
import Taiji.Pipeline.SC.ATACSeq.Functions.QC

plotClusters :: SCATACSeqConfig config
             => SCATACSeq S [CellCluster]
             -> ReaderT config IO ()
plotClusters input = do
    dir <- asks ((<> "/Cluster") . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_cluster.html" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    liftIO $ visualizeCluster output Nothing $ input^.replicates._2.files

{-
plotClusters' :: SCATACSeqConfig config
              => [CellCluster]
              -> ReaderT config IO ()
plotClusters' clusters = do
    dir <- asks ((<> "/Cluster") . _scatacseq_output_dir) >>= getPath
    let output = dir <> "/cluster.html"
    liftIO $ visualizeCluster' output clusters
    -}