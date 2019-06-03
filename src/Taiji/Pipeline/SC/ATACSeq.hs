{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
module Taiji.Pipeline.SC.ATACSeq (builder) where

import           Control.Workflow

import           Taiji.Prelude
import           Taiji.Pipeline.SC.ATACSeq.Functions

builder :: Builder ()
builder = do
    node "Read_Input" 'readInput $
        doc .= "Read ATAC-seq data information from input file."
    node "Download_Data" 'downloadData $
        doc .= "Download data."
    node "Get_Fastq" [| return . getFastq |] $ return ()
    nodePar "Align" 'tagAlign $ do
        nCore .= 8
        doc .= "Read alignment using BWA. The default parameters are: " <>
            "bwa mem -M -k 32."
    nodePar "Filter_Bam" 'filterBamSort $ do
        doc .= "Remove low quality tags using: samtools -F 0x70c -q 30"
    nodePar "QC" 'qualityControl $ return ()
    nodePar "Get_Bins" 'getBins $ return ()
    path [ "Read_Input", "Download_Data", "Get_Fastq", "Align", "Filter_Bam"
         , "QC" ]

    node "Get_Bed" [| \(input, x) -> return $ getSortedBed input ++ 
        (traverse.replicates._2.files %~ (^._1) $ x) |] $ return ()
    nodePar "Make_Count_Matrix" 'mkCellByBinMat $ return ()
    [ "Download_Data", "QC"] ~> "Get_Bed"
    nodePar "TF_IDF" 'lsaClust $ return ()
    path ["Get_Bed", "Get_Bins", "Make_Count_Matrix", "TF_IDF"]

    node "Merge_Count_Matrix_Prep" [| \(x, y) -> return $
        zipExp (x & mapped.replicates._2.files %~ (^._2)) y
        |]$ return ()
    node "Merge_Count_Matrix" 'mergeMatrix $ return ()
    node "LSA_Merged" 'lsaClustMerged $ return ()
    ["Get_Bins", "Make_Count_Matrix"] ~> "Merge_Count_Matrix_Prep"
    path ["Merge_Count_Matrix_Prep", "Merge_Count_Matrix", "LSA_Merged"]

    node "Plot_Cluster_Prep" [| \(x,y) -> return $
        zipExp x $ y & mapped.replicates._2.files %~ (^._3)
        |] $ return ()
    nodePar "Plot_Cluster" 'plotClusters $ return ()
    ["TF_IDF", "QC"] ~> "Plot_Cluster_Prep"
    ["Plot_Cluster_Prep"] ~> "Plot_Cluster"

    node "Make_Bed_Cluster_Prep" [| \(x,y) -> return $ zipExp x y |] $ return ()
    nodePar "Make_Bed_Cluster" 'mkCellClusterBed $ return ()
    nodePar "Subsample_Bed_Cluster" 'subSampleClusterBed $ return ()
    node "Call_Peak_Cluster_Prep" [| return . concatMap split |] $ return ()
    nodePar "Call_Peak_Cluster" 'callPeakCluster $ return ()
    ["Get_Bed", "TF_IDF"] ~> "Make_Bed_Cluster_Prep"
    path ["Make_Bed_Cluster_Prep", "Make_Bed_Cluster",
        "Subsample_Bed_Cluster", "Call_Peak_Cluster_Prep", "Call_Peak_Cluster"]

    nodePar "Make_CutSite_Index" 'mkCutSiteIndex $ return ()
    path ["Get_Bed", "Make_CutSite_Index"]

    node "Get_Open_Region" 'getOpenRegion $ return ()
    node "Find_TFBS_Prep" [| findMotifsPre 5e-5 |] $ return ()
    nodePar "Find_TFBS" 'findMotifs $ return ()
    path ["Get_Bed", "Get_Open_Region", "Find_TFBS_Prep", "Find_TFBS"]

    {-
    -- Snap pipeline
    nodePar "Snap_Pre" 'snapPre $ return ()
    nodePar "Snap_Cluster" [| liftIO . getClusters |] $ return ()
    path ["Get_Bed", "Snap_Pre", "Snap_Cluster"]
    -}

    nodePar "Estimate_Gene_Expr" 'estimateExpr $ return ()
    node "Make_Expr_Table" 'mkExprTable $ return ()
    path ["Call_Peak_Cluster_Prep", "Estimate_Gene_Expr", "Make_Expr_Table"]
