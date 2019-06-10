{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
module Taiji.Pipeline.SC.ATACSeq (builder) where

import           Control.Workflow

import           Taiji.Prelude
import           Taiji.Pipeline.SC.ATACSeq.Functions

builder :: Builder ()
builder = do
    -- The basics
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
    path [ "Read_Input", "Download_Data", "Get_Fastq", "Align", "Filter_Bam"
         , "QC" ]

    node "Get_Bed" [| \(input, x) -> return $ getSortedBed input ++ 
        (traverse.replicates._2.files %~ (^._1) $ x) |] $ return ()
    nodePar "Get_Bins" 'getBins $ return ()
    nodePar "Make_Count_Matrix" 'mkWindowMat $ return ()
    [ "Download_Data", "QC"] ~> "Get_Bed"
    path ["Get_Bed", "Get_Bins", "Make_Count_Matrix"]

    -- LSA
    lsaClust "Cluster/LSA/"
    path ["Make_Count_Matrix", "LSA"]

    {-
    node "Merge_Count_Matrix_Prep" [| \(x, y) -> return $
        zipExp (x & mapped.replicates._2.files %~ (^._2)) y
        |]$ return ()
    node "Merge_Count_Matrix" 'mergeMatrix $ return ()
    node "LSA_Merged" 'performLSAMerged $ return ()
    node "Cluster_LSA_Merged" [| \input -> do
        tmp <- asks _scatacseq_temp_dir
        case input of
            Nothing -> return []
            Just f -> liftIO $ clust True tmp f
        |] $ memory .= 20
    node "Plot_Cluster_Merged" [| \input -> if null input
        then return ()
        else plotClusters' input
        |] $ return ()
    ["Get_Bins", "Make_Count_Matrix"] ~> "Merge_Count_Matrix_Prep"
    path ["Merge_Count_Matrix_Prep", "Merge_Count_Matrix", "LSA_Merged",
        "Cluster_LSA_Merged", "Plot_Cluster_Merged"]

    node "Make_Bed_Cluster_Prep"  [| return . getClusterBarcodes |] $ return ()
    ["Get_Bed", "Cluster_LSA_Merged"] ~> "Make_Bed_Cluster_Prep"
    nodePar "Make_Bed_Cluster" 'getBedCluster $ return ()
    node "Merge_Bed_Cluster" 'mergeBedCluster $ return ()
    nodePar "Call_Peak_Cluster" 'callPeakCluster $ return ()
    node "Merge_Peaks" 'mergePeaks $ return ()
    path ["Make_Bed_Cluster_Prep", "Make_Bed_Cluster", "Merge_Bed_Cluster",
        "Call_Peak_Cluster", "Merge_Peaks"]
        -}

    {-
    -- LDA
    nodePar "LDA" 'performLDA $ return ()
    nodePar "Cluster_LDA" [| \input -> do
        tmp <- asks _scatacseq_temp_dir
        input & replicates.traversed.files %%~ liftIO . clust False tmp
        |] $ return ()
    nodePar "Visualize_LDA_Cluster" [| \x -> do
        dir <- asks ((<> "/Cluster/LDA") . _scatacseq_output_dir) >>= getPath
        liftIO $ plotClusters dir x
        |] $ return ()
    path ["Make_Count_Matrix", "LDA", "Cluster_LDA", "Visualize_LDA_Cluster"]
    -}

    {-
    node "Make_Bed_Cluster_Prep" [| \(x,y) -> return $ zipExp x y |] $ return ()
    nodePar "Make_Bed_Cluster" 'mkCellClusterBed $ return ()
    nodePar "Subsample_Bed_Cluster" 'subSampleClusterBed $ return ()
    node "Call_Peak_Cluster_Prep" [| return . concatMap split |] $ return ()
    nodePar "Call_Peak_Cluster" 'callPeakCluster $ return ()
    ["Get_Bed", "Cluster_LSA"] ~> "Make_Bed_Cluster_Prep"
    path ["Make_Bed_Cluster_Prep", "Make_Bed_Cluster",
        "Subsample_Bed_Cluster", "Call_Peak_Cluster_Prep", "Call_Peak_Cluster"]
        -}

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

    {-
    nodePar "Estimate_Gene_Expr" 'estimateExpr $ return ()
    node "Make_Expr_Table" 'mkExprTable $ return ()
    path ["Call_Peak_Cluster_Prep", "Estimate_Gene_Expr", "Make_Expr_Table"]
    -}
