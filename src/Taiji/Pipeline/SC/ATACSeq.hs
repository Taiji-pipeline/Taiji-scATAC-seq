{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
{-# LANGUAGE LambdaCase #-}
module Taiji.Pipeline.SC.ATACSeq (builder) where

import           Control.Workflow
import qualified Data.Text as T
import qualified Data.ByteString.Char8 as B
import qualified Data.HashMap.Strict as M
import Data.Binary

import           Taiji.Prelude
import           Taiji.Pipeline.SC.ATACSeq.Functions
import Taiji.Pipeline.SC.ATACSeq.Functions.Utils (concatMatrix)

-- | The basic analysis.
basicAnalysis :: Builder ()
basicAnalysis = do
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
        doc .= "Read alignment using BWA. The default parameters are: " <>
            "bwa mem -M -k 32."
    nodePar "Filter_Bam" 'filterBamSort $ do
        doc .= "Remove low quality tags using: samtools -F 0x70c -q 30"
    path ["Align_Prep", "Align", "Filter_Bam"]

    node "Get_Bam" [| \(input, x) -> return $ getBam input ++ x |] $ return ()
    ["Download_Data", "Filter_Bam"] ~> "Get_Bam"

    nodePar "Remove_Duplicates" 'deDuplicates $ return ()
    nodePar "Filter_Cell" 'filterCell $ return ()
    path ["Get_Bam", "Remove_Duplicates", "Filter_Cell"]
    node "Get_Bed" [| \(input, x) -> return $ getSortedBed input ++ x |] $ return ()
    [ "Download_Data", "Filter_Cell"] ~> "Get_Bed"

    node "Get_Genes" [| \_ -> getGeneNames |] $ return ()

-- PreClustering and doublet detection
preClustering :: Builder ()
preClustering = do
    path ["Get_Bed", "Pre_Get_Windows"]
    ["Get_Bed", "Pre_DM_Cluster"] ~> "Pre_Extract_Tags_Prep"
    ["Pre_Make_Peak_Mat", "Remove_Duplicates"] ~> "Pre_Detect_Doublet_Prep"
    ["Pre_Get_Windows", "Get_Genes"] ~> "Pre_Make_Gene_Mat_Prep"
    ["Get_Genes"] ~> "Pre_Get_Markers"
    namespace "Pre" $ do
        -- Creating Cell by Window matrix
        nodePar "Get_Windows" [| getWindows "/temp/Pre/Window/" |] $ return ()
        nodePar "Make_Window_Mat" [| mkWindowMat "/temp/Pre/Window/" |] $ return ()
        path ["Get_Windows", "Make_Window_Mat"]

        -- Clustering in each sample
        --dmClust "/temp/Pre/Cluster/" defClustOpt{_normalization = None, _resolution = Just 0.002}
        dmClust "/temp/Pre/Cluster/" defClustOpt{_normalization = None}
        path ["Make_Window_Mat", "DM_Reduce"]

        -- Extract tags for each cluster
        node "Extract_Tags_Prep"  [| return . uncurry zipExp |] $ return ()
        nodePar "Extract_Tags" [| \input -> input & replicates.traverse.files %%~ ( \(bed, cl) -> do
            let idRep = asDir $ "/temp/Pre/Bed/" <> T.unpack (input^.eid) <>
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
        path ["Extract_Tags_Prep", "Extract_Tags"]

        -- Make cell by peak matrix
        nodePar "Call_Peaks" [| \input -> input & replicates.traverse.files %%~ 
            mapM (findPeaks $ "/temp/Pre/Peak/" <> T.unpack (input^.eid) <> "/") 
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
        node "Get_Markers" [| \fl -> do
            dir <- asks _scatacseq_output_dir >>= getPath
            let output = dir <> "/temp/Pre/marker_genes_idx.txt"
            genes <- liftIO $ map (head . B.split '\t') . B.lines <$> B.readFile fl
            let geneIdx = M.fromList $ zip genes [0 :: Int ..]
            markers <- asks _scatacseq_marker_gene_list >>= \case
                Nothing -> return []
                Just m -> liftIO $
                    map ((\x -> B.pack $ show $ M.lookupDefault undefined x geneIdx) . head . B.split '\t') . B.lines <$> B.readFile m
            liftIO $ B.writeFile output $ B.unlines markers
            return output
            |] $ return ()
        node "Make_Gene_Mat_Prep" [| \(xs, genes) -> return $ zip xs $ repeat genes |] $ return ()
        nodePar "Make_Gene_Mat" [| mkCellByGene "/temp/Pre/Gene/" |] $
            doc .= "Create cell by gene matrix for each sample."
        path ["Make_Gene_Mat_Prep", "Make_Gene_Mat"]

        -- Differetial genes
        nodePar "Get_Ref_Cells" [| liftIO . sampleCells 20 |] $ return ()
        ["DM_Cluster"] ~> "Get_Ref_Cells"
        node "Make_Ref_Gene_Mat" [| \(ref, mat) -> do
            dir <- asks _scatacseq_output_dir >>= getPath
            let output = dir <> "/temp/Pre/merged_ref.mat.gz" 
            mats <- mapM (mkRefMat "/temp/Pre/" False) $ zipExp ref mat
            liftIO $ concatMatrix output $ flip map mats $ \x -> (Nothing, x^.replicates._2.files.location)
            return $ location .~ output $ emptyFile
            |] $ return ()
        ["Get_Ref_Cells", "Make_Gene_Mat"] ~> "Make_Ref_Gene_Mat"
        node "Extract_Cluster_Gene_Matrix" [| \(ref, mat, cl, marker) -> do
            res <- mapM (extractSubMatrix "/temp/Pre/Cluster/") $ zipExp mat cl
            return $ zip3 (concat res) (repeat ref) $ repeat marker
            |] $ return ()
        ["Make_Ref_Gene_Mat", "Make_Gene_Mat", "DM_Cluster", "Get_Markers"] ~> "Extract_Cluster_Gene_Matrix"

        nodePar "Diff_Gene" [| \(x,y,z) -> diffGenes "/temp/Pre/Diff/" (Just z) (x,y) |] $ return ()
        node "Marker_Enrichment" [| computeMarkerEnrichment "/temp/Pre/Diff/" |] $ return ()
        path ["Extract_Cluster_Gene_Matrix", "Diff_Gene", "Marker_Enrichment"]

        -- Doublet detection
        node "Detect_Doublet_Prep" [| return . uncurry zipExp |] $ return ()
        nodePar "Detect_Doublet" 'detectDoublet $ return ()
        path ["Detect_Doublet_Prep", "Detect_Doublet"]

        node "Cluster_QC_Prep" [| \(x1,x2,x3,x4) -> return $
            (zipExp x1 $ zipExp x2 $ zipExp x3 x4) &
                traverse.replicates.traverse.files %~ (\(a,(b,(c,d))) -> (a,b,c,d))
            |] $ return ()
        nodePar "Cluster_QC" 'plotClusterQC $ return ()
        ["Detect_Doublet", "DM_Cluster", "Make_Gene_Mat", "Marker_Enrichment"]
            ~> "Cluster_QC_Prep"
        ["Cluster_QC_Prep"] ~> "Cluster_QC"

builder :: Builder ()
builder = do
    basicAnalysis
    preClustering

    node "Remove_Doublets_Prep" [| return . uncurry zipExp |] $ return ()
    ["Get_Bed", "Pre_Detect_Doublet"] ~> "Remove_Doublets_Prep"
    nodePar "Remove_Doublets" 'removeDoublet $ return ()
    path ["Remove_Doublets_Prep", "Remove_Doublets"]
--------------------------------------------------------------------------------
-- QC
--------------------------------------------------------------------------------
    node "QC" 'plotStat $ return ()
    ["Remove_Duplicates"] ~> "QC"

--------------------------------------------------------------------------------
-- Creating Cell by Window matrix
--------------------------------------------------------------------------------
    nodePar "Get_Windows" [| getWindows "/Feature/Window/" |] $ return ()
    nodePar "Make_Window_Matrix" [| mkWindowMat "/Feature/Window/" |] $ return ()
    path ["Get_Bed", "Get_Windows", "Make_Window_Matrix"]

    -- merged matrix
    node "Merge_Window_Matrix_Prep" [| \(x, y) -> return $
        zipExp (x & mapped.replicates._2.files %~ (^._2)) y
        |]$ return ()
    node "Merge_Window_Matrix" [| mergeFeatMatrix "/Feature/Window/merged_cell_by_window" |] $ return ()
    ["Get_Windows", "Make_Window_Matrix"] ~> "Merge_Window_Matrix_Prep"
    path ["Merge_Window_Matrix_Prep", "Merge_Window_Matrix"]

--------------------------------------------------------------------------------
-- Clustering
--------------------------------------------------------------------------------
    lsaBuilder

    node "Combine_SubCluster" [| fmap return . combineClusters "/Cluster_by_peak/" |] $ return ()
    path ["SubCluster_LSA_Cluster", "Combine_SubCluster"]

--------------------------------------------------------------------------------
-- Creating Cell by Peak matrix
--------------------------------------------------------------------------------
    -- Extract tags for each cluster
    extractTags "/Bed/Cluster/"
    ["Get_Bed", "Peak_LSA_Cluster"] ~> "Extract_Tags_Prep"

    -- Extract tags for subclusters
    namespace "SubCluster" $ extractTags "/Bed/SubCluster/"
    ["Get_Bed", "SubCluster_LSA_Cluster"] ~> "SubCluster_Extract_Tags_Prep"

    -- Call peaks in each cluster and generate peak matrix
    genPeakMat "/Feature/Peak/Cluster/" Nothing "Merge_Tags" "Get_Windows"
    node "Extract_Cluster_Matrix" [| \(x,y) -> 
        let [input] = zipExp x y
        in extractSubMatrix "/Feature/Peak/Cluster/" input |] $ return ()
    ["Merge_Peak_Matrix", "Peak_LSA_Cluster"] ~> "Extract_Cluster_Matrix"

    node "Cluster_Correlation" 'clusterCorrelation $ return ()
    ["Call_Peaks", "Merge_Peaks"] ~> "Cluster_Correlation"

    -- Call peaks SubCluster
    genPeakMat "/Feature/Peak/SubCluster/" (Just "SubCluster")
        "SubCluster_Merge_Tags" "Get_Windows"

    node "SubCluster_Correlation" 'clusterCorrelation $ return ()
    ["SubCluster_Call_Peaks", "SubCluster_Merge_Peaks"] ~> "SubCluster_Correlation"

--------------------------------------------------------------------------------
-- Creating Cell by Gene matrix
--------------------------------------------------------------------------------
    node "Make_Gene_Mat_Prep" [| \(xs, genes) -> return $ zip xs $ repeat genes |] $ return ()
    nodePar "Make_Gene_Mat" [| mkCellByGene "/Feature/Gene/" |] $
        doc .= "Create cell by gene matrix for each sample."
    ["Get_Windows", "Get_Genes"] ~> "Make_Gene_Mat_Prep"
    node "Merge_Gene_Mat" 'mergeCellByGeneMatrix $ return ()
    path ["Make_Gene_Mat_Prep", "Make_Gene_Mat", "Merge_Gene_Mat"]

    node "Extract_Cluster_Gene_Matrix" [| \(x,y) -> 
        let [input] = zipExp x y
        in extractSubMatrix "/Feature/Gene/Cluster/" input |] $ return ()
    ["Merge_Gene_Mat", "Peak_LSA_Cluster"] ~> "Extract_Cluster_Gene_Matrix"

    -- SubClusters
    node "Extract_SubCluster_Gene_Matrix" [| \(x,y) ->
        let [input] = zipExp x y
        in extractSubMatrix "/Feature/Gene/SubCluster/" input |] $ return ()
    ["Merge_Gene_Mat", "Combine_SubCluster"] ~> "Extract_SubCluster_Gene_Matrix"

--------------------------------------------------------------------------------
-- Differential Peak analysis
--------------------------------------------------------------------------------
    node "Get_Ref_Cells" [| \[x] -> liftIO $ sampleCells 200 x |] $ return ()
    ["Peak_LSA_Cluster"] ~> "Get_Ref_Cells"

    node "Get_SubCluster_Ref_Cells" [| \[x] -> liftIO $ sampleCells 100 x |] $ return ()
    ["Combine_SubCluster"] ~> "Get_SubCluster_Ref_Cells"

    node "Make_Ref_Peak_Mat" [| \(ref, mat) -> do
        let [input] = zipExp [ref] mat
        mkRefMat "/Feature/Peak/" False input
        |] $ return ()
    ["Get_Ref_Cells", "Merge_Peak_Matrix"] ~> "Make_Ref_Peak_Mat"

    node "Diff_Peak_Prep" [| \(x, ref) -> return $ zip x $ repeat $ ref^.replicates._2.files |] $ return ()
    ["Extract_Cluster_Matrix", "Make_Ref_Peak_Mat"] ~> "Diff_Peak_Prep"
    nodePar "Diff_Peak" 'diffPeaks $ return ()
    path ["Diff_Peak_Prep", "Diff_Peak"]

    node "RPKM_Peak_Prep" [| \(x, y) -> return $ zip x $ repeat y |] $ return ()
    nodePar "RPKM_Peak" 'rpkmPeak $ return ()
    ["Merge_Tags", "Merge_Peaks"] ~> "RPKM_Peak_Prep"
    path ["RPKM_Peak_Prep", "RPKM_Peak"]
    node "RPKM_Diff_Peak" 'rpkmDiffPeak $ return () 
    ["Diff_Peak", "Merge_Peaks", "RPKM_Peak"] ~> "RPKM_Diff_Peak"

--------------------------------------------------------------------------------
-- Differential Gene analysis
--------------------------------------------------------------------------------
    node "Make_Ref_Gene_Mat" [| \(ref, mat) -> do
        let [input] = zipExp [ref] mat
        mkRefMat "/Feature/Gene/" True input |] $ return ()
    ["Get_Ref_Cells", "Make_Gene_Mat"] ~> "Make_Ref_Gene_Mat"

    node "Diff_Gene_Prep" [| \(x, ref) -> return $ zip x $ repeat $ ref^.replicates._2.files |] $ return ()
    ["Extract_Cluster_Gene_Matrix", "Make_Ref_Gene_Mat"] ~> "Diff_Gene_Prep"
    nodePar "Diff_Gene" [| diffGenes "/Diff/Gene/" Nothing |] $ return ()
    path ["Diff_Gene_Prep", "Diff_Gene"]

    -- SubCluster
    node "Make_SubCluster_Ref_Gene_Mat" [| \(ref, mat) -> do
        let [input] = zipExp [ref] mat
        mkRefMat "/Feature/Gene/Subcluster/" True input |] $ return ()
    ["Get_SubCluster_Ref_Cells", "Make_Gene_Mat"] ~> "Make_SubCluster_Ref_Gene_Mat"

    node "SubCluster_Diff_Gene_Prep" [| \(x, ref) -> return $ zip x $ repeat $ ref^.replicates._2.files |] $ return ()
    ["Extract_SubCluster_Gene_Matrix", "Make_SubCluster_Ref_Gene_Mat"] ~> "SubCluster_Diff_Gene_Prep"
    nodePar "SubCluster_Diff_Gene" [| diffGenes "/Diff/Gene/" Nothing |] $ return ()
    path ["SubCluster_Diff_Gene_Prep", "SubCluster_Diff_Gene"]


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
    -- Snap pipeline
    nodePar "Snap_Pre" 'snapPre $ return ()
    nodePar "Snap_Reduce" 'performSnap $ return ()
    nodePar "Snap_Cluster" [| doClustering "/Snap/" $ defClustOpt{_normalization=None} |] $ return ()
    nodePar "Snap_Viz" [| \x -> do
        dir <- asks ((<> "/Snap/" ) . _scatacseq_output_dir) >>= getPath
        liftIO $ plotClusters dir x
        |] $ return ()
    path ["Get_Bed", "Snap_Pre", "Snap_Reduce", "Snap_Cluster", "Snap_Viz"]

    nodePar "Snap_Mat" 'mkSnapMat $ return ()
    nodePar "Snap_new_Cluster" [| doClustering "/Snap_new/" $ defClustOpt{_normalization=None} |] $ return ()
    nodePar "Snap_new_Viz" [| \x -> do
        dir <- asks ((<> "/Snap_new/" ) . _scatacseq_output_dir) >>= getPath
        liftIO $ plotClusters dir x
        |] $ return ()
    path ["Make_Window_Matrix", "Snap_Mat", "Snap_new_Cluster", "Snap_new_Viz"]

    nodePar "Snap_Merged_Reduce" 'mkSnapMat $ return ()
    nodePar "Snap_Merged_Cluster" [| doClustering "/Cluster_by_peak/Snap/" $ defClustOpt{_normalization=None} |] $ return ()
    nodePar "Snap_Merged_Viz" [| \x -> do
        dir <- asks ((<> "/Cluster_by_peak/Snap/" ) . _scatacseq_output_dir) >>= getPath
        liftIO $ plotClusters dir x
        |] $ return ()
    path ["Merge_Window_Matrix", "Snap_Merged_Reduce", "Snap_Merged_Cluster", "Snap_Merged_Viz"]
    -}