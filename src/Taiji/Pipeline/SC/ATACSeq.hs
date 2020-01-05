{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
{-# LANGUAGE LambdaCase #-}
module Taiji.Pipeline.SC.ATACSeq (builder) where

import           Control.Workflow
import qualified Data.Text as T
import qualified Data.ByteString.Char8 as B
import           Bio.Seq.IO (withGenome, getChrSizes)
import           Bio.Data.Bed (streamBedGzip)
import Data.Binary

import           Taiji.Prelude
import           Taiji.Utils (concatMatrix)
import           Taiji.Pipeline.SC.ATACSeq.Functions
import Taiji.Pipeline.SC.ATACSeq.Types

-- | The basic analysis.
basicAnalysis :: Builder ()
basicAnalysis = do
    node "Read_Input" 'readInput $
        doc .= "Read ATAC-seq data information from input file."
    node "Download_Data" 'download $ doc .= "Download data."
    node "Get_Fastq" [| return . getFastq |] $ return ()
    node "Make_Index" 'mkIndices $ doc .= "Generate the BWA index."
    path ["Read_Input", "Download_Data", "Get_Fastq", "Make_Index"]
 
    node "Align_Prep" [| return . fst |] $ return ()
    ["Get_Fastq", "Make_Index"] ~> "Align_Prep"
    nodePar "Align" 'tagAlign $ do
        nCore .= 8
        doc .= "Read alignment using BWA. The default parameters are: " <>
            "bwa mem -M -k 32."
    nodePar "Filter_Bam" 'filterNameSortBam $ do
        doc .= "Remove low quality tags using: samtools -F 0x70c -q 30"
    path ["Align_Prep", "Align"]

    node "Filter_Bam_Prep" [| \(input, x) -> return $ getBamUnsorted input ++ x |] $ return ()
    ["Download_Data", "Align"] ~> "Filter_Bam_Prep"
    ["Filter_Bam_Prep"] ~> "Filter_Bam"

    node "Get_Bam" [| \(input, x) -> return $ getBam input ++ x |] $ return ()
    ["Download_Data", "Filter_Bam"] ~> "Get_Bam"

    nodePar "Remove_Duplicates" 'deDuplicates $ return ()
    nodePar "Filter_Cell" 'filterCell $ return ()
    path ["Get_Bam", "Remove_Duplicates", "Filter_Cell"]
    node "Get_Bed" [| \(input, x) -> do
        let output = getSortedBed input ++ x
        unless (null output) $ getGenomeIndex >> return ()
        return output
        |] $ return ()
    [ "Download_Data", "Filter_Cell"] ~> "Get_Bed"

-- PreClustering and doublet detection
preClustering :: Builder ()
preClustering = do
    path ["Get_Bed", "Pre_Get_Windows"]
    ["Get_Bed", "Pre_Cluster"] ~> "Pre_Extract_Tags_Prep"
    ["Pre_Make_Peak_Mat", "Remove_Duplicates"] ~> "Pre_Detect_Doublet_Prep"
    ["Get_Bed", "Pre_Detect_Doublet"] ~> "Pre_Remove_Doublets_Prep"
    namespace "Pre" $ do
        -- Creating Cell by Window matrix
        nodePar "Get_Windows" [| getWindows "/Feature/Window/" |] $ return ()
        nodePar "Make_Window_Mat" [| mkWindowMat "/Feature/Window/" |] $ return ()
        path ["Get_Windows", "Make_Window_Mat"]

        -- Clustering in each sample
        nodePar "Cluster" [| \input -> do
            let prefix = "/temp/Pre/Cluster/"
            filterMatrix prefix input >>= spectral prefix Nothing >>=
                mkKNNGraph prefix >>= clustering prefix 1 RBConfiguration 
            |] $ return ()
        path ["Make_Window_Mat", "Cluster"]

        -- Extract tags for each cluster
        node "Extract_Tags_Prep"  [| return . uncurry zipExp |] $ return ()
        nodePar "Extract_Tags" [| \input -> input & replicates.traverse.files %%~ ( \(bed, cl) -> do
            let idRep = asDir $ "/temp/Pre/Bed/" <> T.unpack (input^.eid) <>
                    "_rep" <> show (input^.replicates._1)
            dir <- asks _scatacseq_output_dir >>= getPath . (<> idRep)
            clusters <- liftIO $ decodeFile $ cl^.location
            runResourceT $ runConduit $ streamBedGzip (bed^.location) .|
                splitBedByCluster' dir clusters
            )
            |] $ return ()
        path ["Extract_Tags_Prep", "Extract_Tags"]

        -- Make cell by peak matrix
        nodePar "Call_Peaks" [| \input -> input & replicates.traverse.files %%~ 
            mapM (findPeaks ("/temp/Pre/Peak/" <> T.unpack (input^.eid) <> "/")
                 ( def & mode .~ NoModel (-100) 200
                       & cutoff .~ QValue 0.05
                 ))
            |] $ return ()
        nodePar "Merge_Peaks" [| \input -> input & replicates.traverse.files %%~ 
            mergePeaks ("/temp/Pre/Peak/" <> T.unpack (input^.eid) <> "/")
            |] $ return ()
        path ["Extract_Tags", "Call_Peaks", "Merge_Peaks"]

        node "Get_Peak_List" [| \inputs -> mergePeaks "/temp/Pre/Peak/" $
            flip concatMap inputs $ \input -> input^.replicates._2.files
            |] $ return ()
        path ["Call_Peaks", "Get_Peak_List"]

        node "Make_Peak_Mat_Prep" [| \(x, y) -> return $ flip map (zipExp x y) $ \input ->
            input & replicates._2.files %~ (\((a,_,c), pk) -> (a,fromJust pk,c))
            |] $ return ()
        nodePar "Make_Peak_Mat" [| mkPeakMat "/temp/Pre/Peak/" |] $ return ()
        ["Get_Windows", "Merge_Peaks"] ~> "Make_Peak_Mat_Prep"
        path ["Make_Peak_Mat_Prep", "Make_Peak_Mat"]

        -- Make cell by gene matrix
        node "Get_Promoters" [| \input -> if null input
            then return Nothing
            else Just <$> writePromoters
            |] $ doc .= "Get the list of promoters from the annotation file."
        ["Get_Windows"] ~> "Get_Promoters"
        node "Make_Gene_Mat_Prep" [| \(xs, genes) -> 
            let xs' = map (\x -> x & replicates.traverse.files %~ (\(a,_,c) -> (a,c))) xs
            in return $ zip xs' $ repeat $ fromJust genes |] $ return ()
        nodePar "Make_Gene_Mat" [| mkCellByGene "/temp/Pre/Gene/" |] $
            doc .= "Create cell by transcript matrix for each sample."
        ["Get_Windows", "Get_Promoters"] ~> "Make_Gene_Mat_Prep"
        path ["Make_Gene_Mat_Prep", "Make_Gene_Mat"]

        -- Doublet detection
        node "Detect_Doublet_Prep" [| return . uncurry zipExp |] $ return ()
        nodePar "Detect_Doublet" 'detectDoublet $ return ()
        path ["Detect_Doublet_Prep", "Detect_Doublet"]
        node "Remove_Doublets_Prep" [| return . uncurry zipExp |] $ return ()
        nodePar "Remove_Doublets" 'removeDoublet $ return ()
        path ["Remove_Doublets_Prep", "Remove_Doublets"]

        -- Make feature matrix
        node "Make_Feat_Mat_Prep" [| \(bed, pk) -> return $
            flip map (zip bed $ repeat $ fromJust pk) $ \(x, p) ->
                x & replicates.traverse.files %~ (\(a,b) -> (a,p,b))
            |] $ return ()
        nodePar "Make_Feat_Mat" [| mkPeakMat "/temp/Pre/Feature/" |] $ return ()
        ["Remove_Doublets", "Get_Peak_List"] ~> "Make_Feat_Mat_Prep"
        node "Merge_Feat_Mat" [| \mats -> if null mats
            then return Nothing
            else do
                dir <- asks _scatacseq_output_dir >>= getPath . (<> (asDir "/temp/Pre/Feature/"))
                let output = dir <> "Merged_cell_by_peak.mat.gz"
                liftIO $ concatMatrix output $ flip map mats $ \mat ->
                    ( Just $ B.pack $ T.unpack $ mat^.eid
                    , mat^.replicates._2.files.location )
                return $ Just $ (head mats & eid .~ "Merged") &
                    replicates._2.files.location .~ output
            |] $ return ()
        path ["Make_Feat_Mat_Prep", "Make_Feat_Mat", "Merge_Feat_Mat"]

        {-
        node "Cluster_QC_Prep" [| \(genes, x1, x2, x3, diff) -> do
            let diff' = concatMap split $ mergeExp $ flip map diff $ \x ->
                    let (i, cl) = T.breakOn "+" $ x^.eid
                    in eid .~ i $ replicates._2.files %~ (,) (T.tail cl) $ x
            return $ zip (repeat genes) $ (zipExp x1 $ zipExp x2 $ zipExp x3 diff') &
                traverse.replicates.traverse.files %~ (\(a,(b,(c,d))) -> (a,b,c,d))
            |] $ return ()
        nodePar "Cluster_QC" 'plotClusterQC $ return ()
        ["Get_Genes", "Detect_Doublet", "Cluster", "Make_Gene_Mat", "Diff_Gene"]
            ~> "Cluster_QC_Prep"
        ["Cluster_QC_Prep"] ~> "Cluster_QC"
        -}

builder :: Builder ()
builder = do
    basicAnalysis
    preClustering

--------------------------------------------------------------------------------
-- QC
--------------------------------------------------------------------------------
    node "QC" 'plotStat $ return ()
    ["Pre_Detect_Doublet"] ~> "QC"

--------------------------------------------------------------------------------
-- Clustering
--------------------------------------------------------------------------------
    node "Merged_Filter_Mat" [| \case
        Nothing -> return Nothing
        Just x -> Just <$> filterMatrix "/Cluster/" x
        |] $ return ()
    node "Merged_Reduce_Dims" [| \case
        Nothing -> return Nothing
        Just x -> Just <$>
            spectral "/Cluster/" (Just 3943) x
        |] $ return ()
    node "Merged_Make_KNN" [| \case
        Nothing -> return Nothing
        Just x -> Just <$> mkKNNGraph "/Cluster/" x
        |] $ nCore .= 4
    node "Merged_Cluster" [| \case
        Nothing -> return Nothing
        Just input -> do
            optimizer <- asks _scatacseq_cluster_optimizer
            resolution <- asks _scatacseq_cluster_resolution
            Just <$> clustering "/Cluster/" resolution optimizer input
        |] $ return ()
    path ["Pre_Merge_Feat_Mat", "Merged_Filter_Mat", "Merged_Reduce_Dims",
        "Merged_Make_KNN", "Merged_Cluster"]
    node "Merged_Cluster_Viz" [| \(qc, input) -> case input of
        Nothing -> return ()
        Just x -> do
            dir <- figDir
            liftIO $ plotClusters dir (qc, x)
        |] $ return ()
    ["QC", "Merged_Cluster"] ~> "Merged_Cluster_Viz"

--------------------------------------------------------------------------------
-- Make Cluster BED and BigWig files
--------------------------------------------------------------------------------
    -- Extract tags for each cluster
    node "Extract_Tags_Prep" [| \(x,y) -> return $ zip x $ repeat $ fromJust y |] $ return ()
    extractTags "/Bed/Cluster/"
    ["Pre_Remove_Doublets", "Merged_Cluster"] ~> "Extract_Tags_Prep"
    ["Extract_Tags_Prep"] ~> "Extract_Tags"

    nodePar "Make_BigWig" [| \(nm, fl) -> do
        dir <- asks _scatacseq_output_dir >>= getPath . (<> "/BigWig/Cluster/")
        seqIndex <- getGenomeIndex
        let output = dir <> B.unpack nm <> ".bw"
        liftIO $ do
            chrSize <- withGenome seqIndex $ return . getChrSizes
            bedToBigWig output chrSize fl
        |] $ return ()
    ["Merge_Tags"] ~> "Make_BigWig"

 --------------------------------------------------------------------------------
-- Make cell by peak matrix
--------------------------------------------------------------------------------
    nodePar "Call_Peaks" [| \x -> do
        opts <- getCallPeakOpt
        findPeaks "/Feature/Peak/Cluster/" opts x
        |] $ return ()
    node "Merge_Peaks" [| mergePeaks "/Feature/Peak/" |] $ return ()
    path ["Merge_Tags", "Call_Peaks", "Merge_Peaks"]

    node "Make_Peak_Mat_Prep" [| \(bed, pk) -> return $
        flip map (zip bed $ repeat $ fromJust pk) $ \(x, p) ->
            x & replicates.traverse.files %~ (\(a,b) -> (a,p,b))
        |] $ return ()
    nodePar "Make_Peak_Mat" [| mkPeakMat "/Feature/Peak/" |] $ return ()
    ["Pre_Remove_Doublets", "Merge_Peaks"] ~> "Make_Peak_Mat_Prep"
    path ["Make_Peak_Mat_Prep", "Make_Peak_Mat"]

--------------------------------------------------------------------------------
-- Differential Peak analysis
--------------------------------------------------------------------------------
    node "Cluster_Peak_Mat" [| \(mats, cl) -> case cl of
        Nothing -> return []
        Just x -> subMatrix "/Feature/Peak/Cluster/" mats $ x^.replicates._2.files
        |] $ return ()
    ["Make_Peak_Mat", "Merged_Cluster"] ~> "Cluster_Peak_Mat"

    node "Peak_Acc" [| computePeakRAS "/Feature/Peak/Cluster/" |] $ return ()
    ["Merge_Peaks", "Cluster_Peak_Mat"] ~> "Peak_Acc"

--------------------------------------------------------------------------------
-- Make gene matrix
--------------------------------------------------------------------------------
    node "Cluster_Gene_Mat" [| \(mats, cl) -> case cl of
        Nothing -> return []
        Just x -> subMatrix "/Feature/Gene/Cluster/" mats $ x^.replicates._2.files
        |] $ return ()
    ["Pre_Make_Gene_Mat", "Merged_Cluster"] ~> "Cluster_Gene_Mat"
    node "Gene_Acc" [| mkExprTable "/Feature/Gene/Cluster/" |] $ return ()
    ["Pre_Get_Promoters", "Cluster_Gene_Mat"] ~> "Gene_Acc"

--------------------------------------------------------------------------------
-- Call CRE interactions
--------------------------------------------------------------------------------
    -- Motif finding
    node "Find_TFBS_Prep" [| findMotifsPre 1e-5 |] $ return ()
    nodePar "Find_TFBS" 'findMotifs $ return ()
    path ["Merge_Peaks", "Find_TFBS_Prep", "Find_TFBS"]
