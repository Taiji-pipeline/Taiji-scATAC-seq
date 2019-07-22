{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
module Taiji.Pipeline.SC.ATACSeq (builder) where

import           Control.Workflow
import Data.Binary (decodeFile)
import qualified Data.Text as T
import qualified Data.ByteString.Char8 as B

import           Taiji.Prelude
import           Taiji.Pipeline.SC.ATACSeq.Functions

builder :: Builder ()
builder = do
--------------------------------------------------------------------------------
-- The basics
--------------------------------------------------------------------------------
    node "Read_Input" 'readInput $
        doc .= "Read ATAC-seq data information from input file."
    node "Download_Data" 'downloadData $
        doc .= "Download data."
    node "Get_Fastq" [| return . getFastq |] $ return ()
    node "Make_Index" 'mkIndices $ doc .= "Generate the BWA index."
    path ["Read_Input", "Download_Data", "Get_Fastq", "Make_Index"]
 
    node "Align_Prep" [| return . fst |] $ return ()
    ["Get_Fastq", "Make_Index"] ~> "Align_Prep"
    nodePar "Align" 'tagAlign $ do
        nCore .= 4
        memory .= 10
        doc .= "Read alignment using BWA. The default parameters are: " <>
            "bwa mem -M -k 32."
    nodePar "Filter_Bam" 'filterBamSort $ do
        doc .= "Remove low quality tags using: samtools -F 0x70c -q 30"
    nodePar "Remove_Duplicates" 'deDuplicates $ return ()
    nodePar "Filter_Cell" 'filterCell $ return ()
    path ["Align_Prep", "Align", "Filter_Bam", "Remove_Duplicates", "Filter_Cell"]
    node "Get_Bed" [| \(input, x) -> return $ getSortedBed input ++ x |] $ return ()
    [ "Download_Data", "Filter_Cell"] ~> "Get_Bed"

--------------------------------------------------------------------------------
-- QC
--------------------------------------------------------------------------------
    node "QC" [| mapM_ plotStat |] $ return ()
    ["Remove_Duplicates"] ~> "QC"

--------------------------------------------------------------------------------
-- Creating Cell by Window matrix
--------------------------------------------------------------------------------
    nodePar "Get_Bins" 'getBins $ return ()
    nodePar "Make_Window_Matrix" 'mkWindowMat $ return ()
    path ["Get_Bed", "Get_Bins", "Make_Window_Matrix"]

    -- merged matrix
    node "Merge_Window_Matrix_Prep" [| \(x, y) -> return $
        zipExp (x & mapped.replicates._2.files %~ (^._2)) y
        |]$ return ()
    node "Merge_Window_Matrix" [| mergeFeatMatrix "/Feature/Window/merged_cell_by_window.txt.gz" |] $ return ()
    ["Get_Bins", "Make_Window_Matrix"] ~> "Merge_Window_Matrix_Prep"
    path ["Merge_Window_Matrix_Prep", "Merge_Window_Matrix"]

--------------------------------------------------------------------------------
-- Process each sample
--------------------------------------------------------------------------------
    -- Clustering in each sample
    namespace "Each" $ lsaClust "/Cluster_by_window/LSA/Each/" $
        ClustOpt UnitBall UMAP Nothing
    path ["Make_Window_Matrix", "Each_LSA_Reduce"]

    -- Extract tags for each cluster
    node "Each_Extract_Tags_Prep"  [| return . uncurry zipExp |] $ return ()
    nodePar "Each_Extract_Tags" [| \input -> input & replicates.traverse.files %%~ ( \(bed, cl) -> do
        let idRep = asDir $ "/temp/Bed/Each/" <> T.unpack (input^.eid) <>
                "_rep" <> show (input^.replicates._1)
        dir <- asks _scatacseq_output_dir >>= getPath . (<> idRep)
        clusters <- liftIO $ decodeFile $ cl^.location
        let (nm, bcs) = unzip $ flip map clusters $ \c ->
                (_cluster_name c, map _cell_barcode $ _cluster_member c)
            outputs = map (\x -> dir <> "/" <> B.unpack x <> ".bed") nm
        fls <- liftIO $ extractBedByBarcode outputs bcs bed
        return $ zip nm fls
        )
        |] $ return ()
    ["Get_Bed", "Each_LSA_Cluster"] ~> "Each_Extract_Tags_Prep"
    path ["Each_Extract_Tags_Prep", "Each_Extract_Tags"]

    -- Make cell by peak matrix
    nodePar "Each_Call_Peaks" [| \input -> input & replicates.traverse.files %%~ 
        mapM (findPeaks $ "/temp/Peak/Each/" <> T.unpack (input^.eid) <> "/") 
        |] $ return ()
    nodePar "Each_Merge_Peaks" [| \input -> input & replicates.traverse.files %%~ 
        mergePeaks ("/temp/Peak/Each/" <> T.unpack (input^.eid) <> "/")
        |] $ return ()
    path ["Each_Extract_Tags", "Each_Call_Peaks", "Each_Merge_Peaks"]

    node "Each_Make_Peak_Matrix_Prep" [| \(x, y) -> return $ flip map (zipExp x y) $ \input ->
        input & replicates._2.files %~ (\((a,_,c), pk) -> (a,fromJust pk,c))
        |] $ return ()
    nodePar "Each_Make_Peak_Matrix" [| mkPeakMat "/temp/Peak/Each/" |] $ return ()
    ["Get_Bins", "Each_Merge_Peaks"] ~> "Each_Make_Peak_Matrix_Prep"
    path ["Each_Make_Peak_Matrix_Prep", "Each_Make_Peak_Matrix"]

    nodePar "Detect_Doublet" 'detectDoublet $ return ()
    path ["Each_Make_Peak_Matrix", "Detect_Doublet"]


--------------------------------------------------------------------------------
-- DM
--------------------------------------------------------------------------------
    -- Clustering in each sample
    namespace "Each" $ dmClust "/Cluster_by_window/DM/Each/"
    path ["Make_Window_Matrix", "Each_DM_Reduce"]

--------------------------------------------------------------------------------
-- LSA
--------------------------------------------------------------------------------
    -- Clustering 1st round (By Window)
    namespace "Window" $ lsaClust "/Cluster_by_window/LSA/" $
        ClustOpt UnitBall UMAP Nothing
    path ["Merge_Window_Matrix", "Window_LSA_Reduce"]

    -- Extract tags for each cluster
    namespace "Window_LSA" $ extractTags "/temp/Bed/Cluster/"
    ["Get_Bed", "Window_LSA_Cluster"] ~> "Window_LSA_Extract_Tags_Prep"

    -- Call peaks 1st round
    genPeakMat "/temp/Peak/" (Just "LSA_1st") 
        "Window_LSA_Merge_Tags" "Get_Bins"

    -- Clustering 2nd round
    namespace "Peak" $ lsaClust "/Cluster_by_peak/LSA/" $
        ClustOpt UnitBall UMAP Nothing
    path ["LSA_1st_Merge_Peak_Matrix", "Peak_LSA_Reduce"]

    -- Subclustering
    node "Extract_Sub_Matrix" 'extractSubMatrix $ return ()
    ["LSA_1st_Merge_Peak_Matrix", "Peak_LSA_Cluster"] ~> "Extract_Sub_Matrix"
    namespace "SubCluster" $ lsaClust "/Cluster_by_peak/LSA/SubCluster/" $
        ClustOpt UnitBall UMAP $ Just 5
    path ["Extract_Sub_Matrix", "SubCluster_LSA_Reduce"]


--------------------------------------------------------------------------------
-- Creating Cell by Peak matrix
--------------------------------------------------------------------------------
    -- Extract tags for each cluster
    extractTags "/Bed/Cluster/"
    ["Get_Bed", "Peak_LSA_Cluster"] ~> "Extract_Tags_Prep"

    -- Extract tags for subclusters
    namespace "SubCluster" $ extractTags "/Bed/SubCluster/"
    ["Get_Bed", "SubCluster_LSA_Cluster"] ~> "SubCluster_Extract_Tags_Prep"

    -- Call peaks final round
    genPeakMat "/Feature/Peak/Cluster/" Nothing "Merge_Tags" "Get_Bins"

    node "Cluster_Correlation" 'clusterCorrelation $ return ()
    ["Call_Peaks", "Merge_Peaks"] ~> "Cluster_Correlation"

    -- Call peaks SubCluster
    genPeakMat "/Feature/Peak/SubCluster/" (Just "SubCluster")
        "SubCluster_Merge_Tags" "Get_Bins"

    node "SubCluster_Correlation" 'clusterCorrelation $ return ()
    ["SubCluster_Call_Peaks", "SubCluster_Merge_Peaks"] ~> "SubCluster_Correlation"

--------------------------------------------------------------------------------
-- Call CRE interactions
--------------------------------------------------------------------------------

    node "Cicero" 'cicero $ return ()
    ["SubCluster_Merge_Peaks", "Merge_Peak_Matrix"] ~> "Cicero"


    -- Estimate gene expression
    nodePar "Estimate_Gene_Expr" 'estimateExpr $ return ()
    node "Make_Expr_Table" [| mkExprTable "expression_profile.tsv" |] $ return ()
    path ["Merge_Tags", "Estimate_Gene_Expr", "Make_Expr_Table"]

    -- Motif finding
    node "Find_TFBS_Prep" [| findMotifsPre 1e-5 |] $ return ()
    nodePar "Find_TFBS" 'findMotifs $ return ()
    path ["SubCluster_Merge_Peaks", "Find_TFBS_Prep", "Find_TFBS"]

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

    {-
    nodePar "Make_CutSite_Index" 'mkCutSiteIndex $ return ()
    path ["Get_Bed", "Make_CutSite_Index"]

    node "Get_Open_Region" 'getOpenRegion $ return ()
    -}

    -- Snap pipeline
    nodePar "Snap_Pre" 'snapPre $ return ()
    nodePar "Snap_Reduce" 'performSnap $ return ()
    nodePar "Snap_Cluster" [| doClustering "/Snap/" $ ClustOpt None UMAP Nothing |] $ return ()
    nodePar "Snap_Viz" [| \x -> do
        dir <- asks ((<> "/Snap/" ) . _scatacseq_output_dir) >>= getPath
        liftIO $ plotClusters dir x
        |] $ return ()
    path ["Get_Bed", "Snap_Pre", "Snap_Reduce", "Snap_Cluster", "Snap_Viz"]

    nodePar "Snap_Mat" 'mkSnapMat $ return ()
    nodePar "Snap_new_Cluster" [| doClustering "/Snap_new/" $ ClustOpt None UMAP Nothing |] $ return ()
    nodePar "Snap_new_Viz" [| \x -> do
        dir <- asks ((<> "/Snap_new/" ) . _scatacseq_output_dir) >>= getPath
        liftIO $ plotClusters dir x
        |] $ return ()
    path ["Make_Window_Matrix", "Snap_Mat", "Snap_new_Cluster", "Snap_new_Viz"]

    nodePar "Snap_Merged_Reduce" 'mkSnapMat $ return ()
    nodePar "Snap_Merged_Cluster" [| doClustering "/Cluster_by_peak/Snap/" $ ClustOpt None UMAP Nothing |] $ return ()
    nodePar "Snap_Merged_Viz" [| \x -> do
        dir <- asks ((<> "/Cluster_by_peak/Snap/" ) . _scatacseq_output_dir) >>= getPath
        liftIO $ plotClusters dir x
        |] $ return ()
    path ["LSA_1st_Merge_Peak_Matrix", "Snap_Merged_Reduce", "Snap_Merged_Cluster", "Snap_Merged_Viz"]