{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
module Taiji.Pipeline.SC.ATACSeq (builder) where

import           Control.Lens
import Bio.Data.Experiment
import           Scientific.Workflow

import           Taiji.Pipeline.SC.ATACSeq.Functions

builder :: Builder ()
builder = do
    nodeS "Read_Input" 'readInput $ do
        submitToRemote .= Just False
        note .= "Read ATAC-seq data information from input file."
    nodeS "Download_Data" 'downloadData $ do
        submitToRemote .= Just False
        note .= "Download data."
    node' "Get_Fastq" 'getFastq $ submitToRemote .= Just False
    nodePS 1 "Align" 'tagAlign $ do
        remoteParam .= "--ntasks-per-node=8"  -- slurm
        --remoteParam .= "-pe smp 2"  -- sge
        note .= "Read alignment using BWA. The default parameters are: " <>
            "bwa mem -M -k 32."
    nodePS 1 "Filter_Bam" 'filterBamSort $ do
        note .= "Remove low quality tags using: samtools -F 0x70c -q 30"
    nodePS 1 "QC" 'qualityControl $ return ()
    nodePS 1 "Get_Bins" 'getBins $ return ()
    path [ "Read_Input", "Download_Data", "Get_Fastq", "Align", "Filter_Bam"
         , "QC" ]

    node' "Get_Bed" [| \(input, x) -> getSortedBed input ++ 
        (traverse.replicates._2.files %~ (^._1) $ x) |] $
        submitToRemote .= Just False
    nodePS 1 "Make_Count_Matrix" 'mkCountMatrix $ return ()
    [ "Download_Data", "QC"] ~> "Get_Bed"
    path ["Get_Bed", "Get_Bins", "Make_Count_Matrix"]

    nodePS 1 "Make_CutSite_Index" 'mkCutSiteIndex $ return ()
    path ["Get_Bed", "Make_CutSite_Index"]

    nodeS "Get_Open_Region" 'getOpenRegion $ return ()
    nodeS "Find_TFBS_Prep" [| findMotifsPre 5e-5 |] $ return ()
    nodePS 1 "Find_TFBS" 'findMotifs $ return ()
    path ["Get_Bed", "Get_Open_Region", "Find_TFBS_Prep", "Find_TFBS"]


    -- Snap pipeline
    nodePS 1 "Snap_Pre" 'snapPre $ return ()
    nodePS 1 "Snap_Cluster" 'getClusters $ return ()
    path ["Get_Bed", "Snap_Pre", "Snap_Cluster"]

