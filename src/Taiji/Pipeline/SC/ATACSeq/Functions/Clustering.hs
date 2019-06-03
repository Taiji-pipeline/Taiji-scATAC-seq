{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE DataKinds #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Clustering
    ( module Taiji.Pipeline.SC.ATACSeq.Functions.Clustering.SnapTools
    , module Taiji.Pipeline.SC.ATACSeq.Functions.Clustering.LSA
    , plotClusters
    ) where

import qualified Data.ByteString.Char8 as B
import qualified Data.Text as T
import qualified Data.HashMap.Strict as M
   
import Taiji.Pipeline.SC.ATACSeq.Functions.Clustering.SnapTools
import Taiji.Pipeline.SC.ATACSeq.Functions.Clustering.LSA
import Taiji.Prelude
import Taiji.Pipeline.SC.ATACSeq.Types
import Taiji.Pipeline.SC.ATACSeq.Functions.Utils
import Taiji.Pipeline.SC.ATACSeq.Functions.QC

plotClusters :: SCATACSeqConfig config
             => SCATACSeq S ([CellCluster], File '[] 'Tsv)
             -> ReaderT config IO ()
plotClusters input = do
    dir <- asks ((<> "/Cluster") . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_cluster.html" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
        (clusters, stat) = input^.replicates._2.files
    liftIO $ do
        coverages <- M.fromList . map
            ((\x -> (_cell_barcode x, log $ fromIntegral $ _uniq_reads x)) . decodeStat) .
            B.lines <$> B.readFile (stat^.location)
        visualizeCluster output coverages clusters