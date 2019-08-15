{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
{-# LANGUAGE LambdaCase #-}
module Taiji.Pipeline.SC.ATACSeq (builder) where

import           Control.Workflow
import qualified Data.Text as T
import qualified Data.ByteString.Char8 as B
import Data.Binary

import           Taiji.Prelude
import           Taiji.Pipeline.SC.ATACSeq.Functions
import Taiji.Pipeline.SC.ATACSeq.Functions.Utils (concatMatrix)
import Taiji.Pipeline.SC.ATACSeq.Types

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


-- PreClustering and doublet detection
preClustering :: Builder ()
preClustering = do
    path ["Get_Bed", "Pre_Get_Windows"]
    ["Get_Bed", "Pre_DM_Cluster"] ~> "Pre_Extract_Tags_Prep"
    ["Pre_Make_Peak_Mat", "Remove_Duplicates"] ~> "Pre_Detect_Doublet_Prep"
    ["Get_Bed", "Pre_Detect_Doublet"] ~> "Pre_Remove_Doublets_Prep"
    namespace "Pre" $ do
        -- Creating Cell by Window matrix
        nodePar "Get_Windows" [| getWindows "/Feature/Window/" |] $ return ()
        nodePar "Make_Window_Mat" [| mkWindowMat "/Feature/Window/" |] $ return ()
        path ["Get_Windows", "Make_Window_Mat"]

        -- Clustering in each sample
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

        node "Make_Peak_Mat_Prep" [| \(x, y) -> return $ flip map (zipExp x y) $ \input ->
            input & replicates._2.files %~ (\((a,_,c), pk) -> (a,fromJust pk,c))
            |] $ return ()
        nodePar "Make_Peak_Mat" [| mkPeakMat "/temp/Pre/Peak/" |] $ return ()
        ["Get_Windows", "Merge_Peaks"] ~> "Make_Peak_Mat_Prep"
        path ["Make_Peak_Mat_Prep", "Make_Peak_Mat"]

        -- Make cell by gene matrix
        node "Get_Genes" [| \_ -> getGeneNames |] $ return ()
        node "Make_Gene_Mat_Prep" [| \(xs, genes) -> 
            let xs' = map (\x -> x & replicates.traverse.files %~ (\(a,_,c) -> (a,c))) xs
            in return $ zip xs' $ repeat genes |] $ return ()
        nodePar "Make_Gene_Mat" [| mkCellByGene "/temp/Pre/Gene/" |] $
            doc .= "Create cell by gene matrix for each sample."
        ["Get_Windows", "Get_Genes"] ~> "Make_Gene_Mat_Prep"
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
        node "Extract_Cluster_Gene_Matrix" [| \(mat, cl) -> fmap concat $
            mapM (extractSubMatrix "/temp/Pre/Cluster/") $ zipExp mat cl
            |] $ return ()
        ["Make_Gene_Mat", "DM_Cluster"] ~> "Extract_Cluster_Gene_Matrix"

        node "Diff_Gene_Prep" [| \(genes, input, ref) -> return $
            zip3 (repeat genes) input $ repeat ref
            |] $ return ()
        nodePar "Diff_Gene" [| diffGenes "/temp/Pre/Diff/" Nothing |] $ return ()
        ["Get_Genes", "Extract_Cluster_Gene_Matrix", "Make_Ref_Gene_Mat"] ~> "Diff_Gene_Prep"
        path ["Diff_Gene_Prep", "Diff_Gene"]

        -- Doublet detection
        node "Detect_Doublet_Prep" [| return . uncurry zipExp |] $ return ()
        nodePar "Detect_Doublet" 'detectDoublet $ return ()
        path ["Detect_Doublet_Prep", "Detect_Doublet"]
        node "Remove_Doublets_Prep" [| return . uncurry zipExp |] $ return ()
        nodePar "Remove_Doublets" 'removeDoublet $ return ()
        path ["Remove_Doublets_Prep", "Remove_Doublets"]

        node "Cluster_QC_Prep" [| \(genes, x1, x2, x3, diff) -> do
            let diff' = concatMap split $ mergeExp $ flip map diff $ \x ->
                    let (i, cl) = T.breakOn "+" $ x^.eid
                    in eid .~ i $ replicates._2.files %~ (,) (T.tail cl) $ x
            return $ zip (repeat genes) $ (zipExp x1 $ zipExp x2 $ zipExp x3 diff') &
                traverse.replicates.traverse.files %~ (\(a,(b,(c,d))) -> (a,b,c,d))
            |] $ return ()
        nodePar "Cluster_QC" 'plotClusterQC $ return ()
        ["Get_Genes", "Detect_Doublet", "DM_Cluster", "Make_Gene_Mat", "Diff_Gene"]
            ~> "Cluster_QC_Prep"
        ["Cluster_QC_Prep"] ~> "Cluster_QC"

builder :: Builder ()
builder = do
    basicAnalysis
    preClustering

--------------------------------------------------------------------------------
-- QC
--------------------------------------------------------------------------------
    node "QC" 'plotStat $ return ()
    ["Remove_Duplicates"] ~> "QC"

--------------------------------------------------------------------------------
-- Make cell by peak matrix
--------------------------------------------------------------------------------
    node "Get_Peak_List" [| \inputs -> mergePeaks "/Feature/Peak/" $
        flip concatMap inputs $ \input -> input^.replicates._2.files
        |] $ return ()
    path ["Pre_Call_Peaks", "Get_Peak_List"]
    node "Make_Peak_Mat_Prep" [| \(bed, pk) -> return $
        flip map (zip bed $ repeat $ fromJust pk) $ \(x, p) ->
            x & replicates.traverse.files %~ (\(a,b) -> (a,p,b))
        |] $ return ()
    nodePar "Make_Peak_Mat" [| mkPeakMat "/Feature/Peak/" |] $ return ()
    ["Pre_Remove_Doublets", "Get_Peak_List"] ~> "Make_Peak_Mat_Prep"
    node "Merge_Peak_Mat" [| \mats -> do
        dir <- asks _scatacseq_output_dir >>= getPath . (<> (asDir "/Feature/Peak/"))
        let output = dir <> "Merged_cell_by_peak.mat.gz"
        liftIO $ concatMatrix output $ flip map mats $ \mat ->
            ( Just $ B.pack $ T.unpack $ mat^.eid
            , mat^.replicates._2.files.location )
        return $ (head mats & eid .~ "Merged") &
            replicates._2.files.location .~ output
        |] $ return ()
    path ["Make_Peak_Mat_Prep", "Make_Peak_Mat", "Merge_Peak_Mat"]

--------------------------------------------------------------------------------
-- Clustering
--------------------------------------------------------------------------------
    node "Merged_Cluster_Reduce" [| performDM "/Cluster/" |] $ return ()
    node "Merged_Cluster" [| doClustering "/Cluster/" defClustOpt{_normalization = None} |] $ return ()
    node "Merged_Cluster_Viz" [| \x -> do
        dir <- figDir
        liftIO $ plotClusters dir x
        |] $ return ()
    path ["Merge_Peak_Mat", "Merged_Cluster_Reduce", "Merged_Cluster", "Merged_Cluster_Viz"]

    -- Subclustering
    node "Extract_Sub_Matrix" [| \(x,y) -> 
        let [input] = zipExp [x] [y]
        in extractSubMatrix "/Feature/Peak/Cluster/" input |] $ return ()
    ["Merge_Peak_Mat", "Merged_Cluster"] ~> "Extract_Sub_Matrix"
    nodePar "Merged_Subcluster_Reduce" [| performDM "/Subcluster/" |] $ return ()
    nodePar "Merged_Subcluster" [| doClustering "/Subcluster/"
        defClustOpt{_normalization = None, _resolution=Just 0.6} |] $ return ()
    nodePar "Merged_Subcluster_Viz" [| \x -> do
        dir <- figDir
        liftIO $ plotClusters dir x
        |] $ return ()
    path ["Extract_Sub_Matrix", "Merged_Subcluster_Reduce", "Merged_Subcluster", "Merged_Subcluster_Viz"]

    -- Extract tags for each cluster
    node "Extract_Tags_Prep" [| \(x,y) -> return $ zip x $ repeat [y] |] $ return ()
    extractTags "/Bed/Cluster/"
    ["Pre_Remove_Doublets", "Merged_Cluster"] ~> "Extract_Tags_Prep"
    ["Extract_Tags_Prep"] ~> "Extract_Tags"
    
    -- Extract tags for subclusters
    node "Subcluster_Extract_Tags_Prep" [| \(x,y) -> return $ zip x $ repeat y |] $ return ()
    namespace "Subcluster" $ extractTags "/Bed/Subcluster/"
    ["Pre_Remove_Doublets", "Merged_Subcluster"] ~> "Subcluster_Extract_Tags_Prep"
    ["Subcluster_Extract_Tags_Prep"] ~> "Subcluster_Extract_Tags"
    

--------------------------------------------------------------------------------
-- Differential genes
--------------------------------------------------------------------------------
    node "Merge_Gene_Mat" 'mergeCellByGeneMatrix $ return ()
    path ["Pre_Make_Gene_Mat", "Merge_Gene_Mat"]

    node "Get_Ref_Cells" [| liftIO . sampleCells 200 |] $ return ()
    ["Merged_Cluster"] ~> "Get_Ref_Cells"
    node "Make_Ref_Gene_Mat" [| \(ref, mat) ->
        let [input] = zipExp [ref] [mat]
        in mkRefMat "/Feature/Gene/" False input
        |] $ return ()
    ["Get_Ref_Cells", "Merge_Gene_Mat"] ~> "Make_Ref_Gene_Mat"

    node "Extract_Cluster_Gene_Matrix" [| \(x,y) -> 
        extractSubMatrix "/Feature/Gene/Cluster/" $
            y & replicates._2.files %~ (,) (x^.replicates._2.files)
        |] $ return ()
    ["Merge_Gene_Mat", "Merged_Cluster"] ~> "Extract_Cluster_Gene_Matrix"
    node "Diff_Gene_Prep" [| \(genes, input, ref) -> return $
        zip3 (repeat genes) input $ repeat $ ref^.replicates._2.files
        |] $ return ()
    nodePar "Diff_Gene" [| diffGenes "/Diff/Gene/" Nothing |] $ return ()
    node "Diff_Gene_Viz" [| plotDiffGene "diff_gene.html" |] $ return ()
    ["Pre_Get_Genes", "Extract_Cluster_Gene_Matrix", "Make_Ref_Gene_Mat"] ~> "Diff_Gene_Prep"
    path ["Diff_Gene_Prep", "Diff_Gene", "Diff_Gene_Viz"]

    node "Extract_Subcluster_Gene_Matrix" [| \(x, ys) -> 
        mapM (extractSubMatrix "/Feature/Gene/Subcluster/") $ flip map ys $ \y ->
            y & replicates._2.files %~ (,) (x^.replicates._2.files)
        |] $ return ()
    ["Merge_Gene_Mat", "Merged_Subcluster"] ~> "Extract_Subcluster_Gene_Matrix"
    node "Subcluster_Diff_Gene_Prep" [| \(genes, input, ref) -> return $
        zip3 (repeat genes) (concat input) $ repeat $ ref^.replicates._2.files
        |] $ return ()
    nodePar "Subcluster_Diff_Gene" [| diffGenes "/Diff/Gene/Subcluster/" Nothing |] $ return ()
    node "Subcluster_Diff_Gene_Viz" [| plotDiffGene "diff_gene_subcluster.html" |] $ return ()
    ["Pre_Get_Genes", "Extract_Subcluster_Gene_Matrix", "Make_Ref_Gene_Mat"] ~> "Subcluster_Diff_Gene_Prep"
    path ["Subcluster_Diff_Gene_Prep", "Subcluster_Diff_Gene", "Subcluster_Diff_Gene_Viz"]

--------------------------------------------------------------------------------
-- Differential Peak analysis
--------------------------------------------------------------------------------
    node "Make_Ref_Peak_Mat" [| \(ref, mat) -> do
        let [input] = zipExp [ref] [mat]
        mkRefMat "/Feature/Peak/" False input
        |] $ return ()
    ["Get_Ref_Cells", "Merge_Peak_Mat"] ~> "Make_Ref_Peak_Mat"

    node "Diff_Peak_Prep" [| \(pk, x, ref) -> return $
        zip3 (repeat $ fromJust pk) x $ repeat $ ref^.replicates._2.files
        |] $ return ()
    ["Get_Peak_List", "Extract_Sub_Matrix", "Make_Ref_Peak_Mat"] ~> "Diff_Peak_Prep"
    nodePar "Diff_Peak" 'diffPeaks $ return ()
    path ["Diff_Peak_Prep", "Diff_Peak"]

    {-
    node "RPKM_Peak_Prep" [| \(x, y) -> return $ zip x $ repeat y |] $ return ()
    nodePar "RPKM_Peak" 'rpkmPeak $ return ()
    ["Merge_Tags", "Merge_Peaks"] ~> "RPKM_Peak_Prep"
    path ["RPKM_Peak_Prep", "RPKM_Peak"]
    node "RPKM_Diff_Peak" 'rpkmDiffPeak $ return () 
    ["Diff_Peak", "Merge_Peaks", "RPKM_Peak"] ~> "RPKM_Diff_Peak"
    -}

--------------------------------------------------------------------------------
-- Call CRE interactions
--------------------------------------------------------------------------------

    node "Cicero" 'cicero $ return ()
    ["Get_Peak_List", "Merge_Peak_Mat"] ~> "Cicero"

    -- Estimate gene expression
    nodePar "Estimate_Gene_Expr" 'estimateExpr $ return ()
    node "Make_Expr_Table" [| mkExprTable "expression_profile.tsv" |] $ return ()
    path ["Merge_Tags", "Estimate_Gene_Expr", "Make_Expr_Table"]

    -- Motif finding
    node "Find_TFBS_Prep" [| findMotifsPre 1e-5 |] $ return ()
    nodePar "Find_TFBS" 'findMotifs $ return ()
    path ["Get_Peak_List", "Find_TFBS_Prep", "Find_TFBS"]