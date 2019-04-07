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
    path ["Read_Input", "Download_Data", "Get_Fastq"]

    {-
    node' "Align_Prep" [| fst |] $ submitToRemote .= Just False
    nodePS 1 "Align" 'atacAlign $ do
        remoteParam .= "--ntasks-per-node=2"  -- slurm
        --remoteParam .= "-pe smp 2"  -- sge
        note .= "Read alignment using BWA. The default parameters are: " <>
            "bwa mem -M -k 32."
    nodePS 1 "Filter_Bam" 'atacFilterBamSort $ do
        note .= "Remove low quality tags using: samtools -F 0x70c -q 30"
    path ["Read_Input", "Download_Data", "Get_Fastq", "Make_Index"]
    ["Get_Fastq", "Make_Index"] ~> "Align_Prep"
    path ["Align_Prep", "Align", "Filter_Bam", "Remove_Duplicates"]
    nodePS 1 "Remove_Duplicates" 'scAtacDeDup $ return ()

    -}

