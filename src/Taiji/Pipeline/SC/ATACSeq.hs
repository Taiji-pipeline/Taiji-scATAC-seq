{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
{-# LANGUAGE LambdaCase #-}
module Taiji.Pipeline.SC.ATACSeq (builder) where

import           Control.Workflow
import qualified Data.Text as T
import qualified Data.ByteString.Char8 as B
import Data.Binary
import System.Random.MWC (uniformR, create)

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
    ["Get_Bed", "Pre_Cluster"] ~> "Pre_Extract_Tags_Prep"
    ["Pre_Make_Peak_Mat", "Remove_Duplicates"] ~> "Pre_Detect_Doublet_Prep"
    ["Get_Bed", "Pre_Detect_Doublet"] ~> "Pre_Remove_Doublets_Prep"
    namespace "Pre" $ do
        -- Creating Cell by Window matrix
        nodePar "Get_Windows" [| getWindows "/Feature/Window/" |] $ return ()
        nodePar "Make_Window_Mat" [| mkWindowMat "/Feature/Window/" |] $ return ()
        path ["Get_Windows", "Make_Window_Mat"]

        -- Clustering in each sample
        spectralClust "/temp/Pre/Cluster/" defClustOpt
        path ["Make_Window_Mat", "Filter_Mat"]

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
        ["Cluster"] ~> "Get_Ref_Cells"
        node "Make_Ref_Gene_Mat" [| \(refs, mats) ->
            let bcs = flip concatMap refs $ \ref ->
                    map (\x -> B.pack (T.unpack $ ref^.eid) <> "+" <> x) $
                        concat $ ref^.replicates._2.files
            in mkRefMat "/temp/Pre/" bcs mats
            |] $ return ()
        ["Get_Ref_Cells", "Make_Gene_Mat"] ~> "Make_Ref_Gene_Mat"
        node "Extract_Cluster_Gene_Matrix" [| \(mat, cl) -> fmap concat $
            mapM (extractSubMatrix "/temp/Pre/Cluster/") $ zipExp mat cl
            |] $ return ()
        ["Make_Gene_Mat", "Cluster"] ~> "Extract_Cluster_Gene_Matrix"

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

        -- Make feature matrix
        node "Make_Feat_Mat_Prep" [| \(bed, pk) -> return $
            flip map (zip bed $ repeat $ fromJust pk) $ \(x, p) ->
                x & replicates.traverse.files %~ (\(a,b) -> (a,p,b))
            |] $ return ()
        nodePar "Make_Feat_Mat" [| mkPeakMat "/temp/Pre/Feature/" |] $ return ()
        ["Remove_Doublets", "Get_Peak_List"] ~> "Make_Feat_Mat_Prep"
        node "Merge_Feat_Mat" [| \mats -> do
            dir <- asks _scatacseq_output_dir >>= getPath . (<> (asDir "/temp/Pre/Feature/"))
            let output = dir <> "Merged_cell_by_peak.mat.gz"
            liftIO $ concatMatrix output $ flip map mats $ \mat ->
                ( Just $ B.pack $ T.unpack $ mat^.eid
                , mat^.replicates._2.files.location )
            return $ (head mats & eid .~ "Merged") &
                replicates._2.files.location .~ output
            |] $ return ()
        path ["Make_Feat_Mat_Prep", "Make_Feat_Mat", "Merge_Feat_Mat"]

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
    node "Merged_Filter_Mat" [| filterMatrix "/Cluster/" |] $ return ()
    node "Merged_Reduce_Dims_Prep" [| \x -> liftIO $ do
        g <- create
        seeds <- replicateM 5 $ uniformR (1::Int, 100000) g
        return $ zip seeds $ repeat x
        |] $ return ()
    nodePar "Merged_Reduce_Dims" [| \(s,x) ->
        spectral ("/Cluster/" ++ show s ++ "/") (Just s) x
        |] $ return ()
    node "Merged_Cluster" [| \x ->
        let [x'] = concatMap split $ mergeExp x
        in clustering "/Cluster/" defClustOpt x'
        |] $ return ()
    path ["Pre_Merge_Feat_Mat", "Merged_Filter_Mat", "Merged_Reduce_Dims_Prep", "Merged_Reduce_Dims",
        "Merged_Cluster"]
    node "Merged_Cluster_Viz" [| \x -> do
        dir <- figDir
        liftIO $ plotClusters dir x
        |] $ return ()
    ["QC", "Merged_Cluster"] ~> "Merged_Cluster_Viz"

    -- Subclustering
    node "Extract_Sub_Matrix" [| \(x,y) -> 
        let [input] = zipExp [x] [y]
        in extractSubMatrix "/temp/Feature/Cluster/" input |] $ return ()
    ["Pre_Merge_Feat_Mat", "Merged_Cluster"] ~> "Extract_Sub_Matrix"

    namespace "Merged_Iterative" $
        spectralClust "/Subcluster/" defClustOpt{_resolution=Just 0.5}
    path ["Extract_Sub_Matrix", "Merged_Iterative_Filter_Mat"]
    node "Merged_Iterative_Cluster_Viz" [| \(qc, xs) -> do
        dir <- figDir
        liftIO $ mapM_ (\x -> plotClusters dir (qc, x)) xs
        |] $ return ()
    ["QC", "Merged_Iterative_Cluster"] ~> "Merged_Iterative_Cluster_Viz"

    -- combine subclusters into a single file
    node "Combine_Clusters" [| \inputs -> do
        dir <- asks _scatacseq_output_dir >>= getPath . (<> "/Subcluster/")
        let output = dir <> "all_subclusters.bin"
        clusters <- fmap concat $ forM inputs $ \input -> liftIO $ do
            cls <- decodeFile $ input^.replicates._2.files.location
            let nm = B.pack $ T.unpack $ input^.eid
            return $ case cls of
                [cl] -> [cl{_cluster_name = nm}]
                xs -> map (\x -> x{_cluster_name = nm <> "." <> B.drop 1 (_cluster_name x)}) xs
        liftIO $ encodeFile output clusters
        return $ head inputs & eid .~ "Subcluster" & replicates._2.files.location .~ output
        |] $ return ()
    ["Merged_Iterative_Cluster"] ~> "Combine_Clusters"

--------------------------------------------------------------------------------
-- Make Cluster BED file
--------------------------------------------------------------------------------
    -- Extract tags for each cluster
    node "Extract_Tags_Prep" [| \(x,y) -> return $ zip x $ repeat y |] $ return ()
    extractTags "/Bed/Cluster/"
    ["Pre_Remove_Doublets", "Merged_Cluster"] ~> "Extract_Tags_Prep"
    ["Extract_Tags_Prep"] ~> "Extract_Tags"
    
    -- Extract tags for subclusters
    node "Subcluster_Extract_Tags_Prep" [| \(x,y) -> return $ zip x $ repeat y |] $ return ()
    namespace "Subcluster" $ extractTags "/Bed/Subcluster/"
    ["Pre_Remove_Doublets", "Combine_Clusters"] ~> "Subcluster_Extract_Tags_Prep"
    ["Subcluster_Extract_Tags_Prep"] ~> "Subcluster_Extract_Tags"

 --------------------------------------------------------------------------------
-- Make cell by peak matrix
--------------------------------------------------------------------------------
    nodePar "Call_Peaks" [| findPeaks "/Feature/Peak/" |] $ return ()
    node "Merge_Peaks" [| mergePeaks "/Feature/Peak/" |] $ return ()
    path ["Subcluster_Merge_Tags", "Call_Peaks", "Merge_Peaks"]

    node "Make_Peak_Mat_Prep" [| \(bed, pk) -> return $
        flip map (zip bed $ repeat $ fromJust pk) $ \(x, p) ->
            x & replicates.traverse.files %~ (\(a,b) -> (a,p,b))
        |] $ return ()
    nodePar "Make_Peak_Mat" [| mkPeakMat "/Feature/Peak/" |] $ return ()
    ["Pre_Remove_Doublets", "Merge_Peaks"] ~> "Make_Peak_Mat_Prep"
    path ["Make_Peak_Mat_Prep", "Make_Peak_Mat"]


    node "Get_Ref_Cells" [| liftIO . sampleCells 200 |] $ return ()
    ["Merged_Cluster"] ~> "Get_Ref_Cells"

--------------------------------------------------------------------------------
-- Differential Peak analysis
--------------------------------------------------------------------------------
    node "Make_Ref_Peak_Mat" [| \(ref, mats) ->
        mkRefMat "/Feature/Peak/" (concat $ ref^.replicates._2.files) mats
        |] $ return ()
    ["Get_Ref_Cells", "Make_Peak_Mat"] ~> "Make_Ref_Peak_Mat"

    {-
    node "Subcluster_Peak_Mat" [| \(mats, cls) -> do
        dir <- asks ((<> "/Feature/Peak/Subcluster/" ) . _scatacseq_output_dir) >>= getPath



        in extractSubMatrix "/Feature/Peak/Subcluster/" input |] $ return ()
    ["Make_Peak_Mat", "Merged_Iterative_Cluster"] ~> "Subcluster_Peak_Mat"
    -}

{-
    node "Diff_Peak_Prep" [| \(pk, x, ref) -> return $
        zip3 (repeat $ fromJust pk) x $ repeat $ ref^.replicates._2.files
        |] $ return ()
    ["Get_Peak_List", "Extract_Sub_Matrix", "Make_Ref_Peak_Mat"] ~> "Diff_Peak_Prep"
    nodePar "Diff_Peak" 'diffPeaks $ return ()
    path ["Diff_Peak_Prep", "Diff_Peak"]
    -}

--------------------------------------------------------------------------------
-- Differential genes
--------------------------------------------------------------------------------
    node "Make_Ref_Gene_Mat" [| \(ref, mat) ->
        mkRefMat "/Feature/Gene/" (concat $ ref^.replicates._2.files) mat
        |] $ return ()
    ["Get_Ref_Cells", "Pre_Make_Gene_Mat"] ~> "Make_Ref_Gene_Mat"

    node "Merge_Gene_Mat" 'mergeCellByGeneMatrix $ return ()
    path ["Pre_Make_Gene_Mat", "Merge_Gene_Mat"]

    node "Extract_Cluster_Gene_Matrix" [| \(x,y) -> 
        extractSubMatrix "/Feature/Gene/Cluster/" $
            y & replicates._2.files %~ (,) (x^.replicates._2.files)
        |] $ return ()
    ["Merge_Gene_Mat", "Merged_Cluster"] ~> "Extract_Cluster_Gene_Matrix"
    node "Diff_Gene_Prep" [| \(genes, input, ref) -> return $
        zip3 (repeat genes) input $ repeat ref
        |] $ return ()
    nodePar "Diff_Gene" [| diffGenes "/Diff/Gene/" Nothing |] $ return ()
    node "Diff_Gene_Viz" [| plotDiffGene "diff_gene.html" |] $ return ()
    ["Pre_Get_Genes", "Extract_Cluster_Gene_Matrix", "Make_Ref_Gene_Mat"] ~> "Diff_Gene_Prep"
    path ["Diff_Gene_Prep", "Diff_Gene", "Diff_Gene_Viz"]

    node "Extract_Subcluster_Gene_Matrix" [| \(x, ys) -> 
        mapM (extractSubMatrix "/Feature/Gene/Subcluster/") $ flip map ys $ \y ->
            y & replicates._2.files %~ (,) (x^.replicates._2.files)
        |] $ return ()
    ["Merge_Gene_Mat", "Merged_Iterative_Cluster"] ~> "Extract_Subcluster_Gene_Matrix"
    node "Subcluster_Diff_Gene_Prep" [| \(genes, input, ref) -> return $
        zip3 (repeat genes) (concat input) $ repeat ref
        |] $ return ()
    nodePar "Subcluster_Diff_Gene" [| diffGenes "/Diff/Gene/Subcluster/" Nothing |] $ return ()
    node "Subcluster_Diff_Gene_Viz" [| plotDiffGene "diff_gene_subcluster.html" |] $ return ()
    ["Pre_Get_Genes", "Extract_Subcluster_Gene_Matrix", "Make_Ref_Gene_Mat"] ~> "Subcluster_Diff_Gene_Prep"
    path ["Subcluster_Diff_Gene_Prep", "Subcluster_Diff_Gene", "Subcluster_Diff_Gene_Viz"]

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

{-
    node "Cicero" 'cicero $ return ()
    ["Get_Peak_List", "Merge_Peak_Mat"] ~> "Cicero"
    -}

    -- Estimate gene expression
    nodePar "Estimate_Gene_Expr" 'estimateExpr $ return ()
    node "Make_Expr_Table" [| mkExprTable "expression_profile.tsv" |] $ return ()
    path ["Merge_Tags", "Estimate_Gene_Expr", "Make_Expr_Table"]

    -- Motif finding
    {-
    node "Find_TFBS_Prep" [| findMotifsPre 1e-5 |] $ return ()
    nodePar "Find_TFBS" 'findMotifs $ return ()
    path ["Get_Peak_List", "Find_TFBS_Prep", "Find_TFBS"]
    -}