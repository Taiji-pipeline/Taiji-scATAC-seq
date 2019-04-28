{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
module Taiji.Pipeline.SC.ATACSeq (builder) where

import           Control.Lens
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
    nodePS 1 "Remove_Duplicates" 'rmDuplicates $ return ()
    nodePS 1 "Count_Tags" 'countTagsMerged $ return ()
    path [ "Read_Input", "Download_Data", "Get_Fastq", "Align", "Filter_Bam"
         , "Remove_Duplicates", "Count_Tags" ]

    nodePS 1 "Make_CutSite_Index" 'mkCutSiteIndex $ return ()
    path ["Remove_Duplicates", "Make_CutSite_Index"]