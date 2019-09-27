{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
{-# LANGUAGE LambdaCase #-}
module Taiji.Pipeline.SC.ATACSeq (builder) where

import           Control.Workflow
import qualified Data.Text as T
import qualified Data.ByteString.Char8 as B
import           Bio.Seq.IO (withGenome, getChrSizes)
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
    node "Download_Data" 'download $ doc .= "Download data."
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
        node "Get_Promoters" [| \_ -> writePromoters |] $
            doc .= "Get the list of promoters from the annotation file."
        node "Make_Transcript_Mat_Prep" [| \(xs, genes) -> 
            let xs' = map (\x -> x & replicates.traverse.files %~ (\(a,_,c) -> (a,c))) xs
            in return $ zip xs' $ repeat genes |] $ return ()
        nodePar "Make_Transcript_Mat" [| mkCellByGene "/temp/Pre/Gene/" |] $
            doc .= "Create cell by transcript matrix for each sample."
        ["Get_Windows", "Get_Promoters"] ~> "Make_Transcript_Mat_Prep"
        path ["Make_Transcript_Mat_Prep", "Make_Transcript_Mat"]

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
    node "Extract_Sub_Matrix" [| \(mats, cls) -> 
        subMatrix "/temp/Feature/Cluster/" mats $ cls^.replicates._2.files
        |] $ return ()
    ["Pre_Make_Feat_Mat", "Merged_Cluster"] ~> "Extract_Sub_Matrix"
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
-- Make Cluster BED and BigWig files
--------------------------------------------------------------------------------
    -- Extract tags for each cluster
    node "Extract_Tags_Prep" [| \(x,y) -> return $ zip x $ repeat y |] $ return ()
    extractTags "/Bed/Cluster/"
    ["Pre_Remove_Doublets", "Merged_Cluster"] ~> "Extract_Tags_Prep"
    ["Extract_Tags_Prep"] ~> "Extract_Tags"

    nodePar "Call_Peaks_Cluster" [| \x -> do
        opts <- asks _scatacseq_callpeak_opts
        findPeaks "/Feature/Peak/Cluster/" opts x
        |] $ return ()
    path ["Merge_Tags", "Call_Peaks_Cluster"]

    -- Extract tags for subclusters
    node "Subcluster_Extract_Tags_Prep" [| \(x,y) -> return $ zip x $ repeat y |] $ return ()
    namespace "Subcluster" $ extractTags "/Bed/Subcluster/"
    ["Pre_Remove_Doublets", "Combine_Clusters"] ~> "Subcluster_Extract_Tags_Prep"
    ["Subcluster_Extract_Tags_Prep"] ~> "Subcluster_Extract_Tags"

    nodePar "Subcluster_Make_BigWig" [| \(nm, fl) -> do
        dir <- asks _scatacseq_output_dir >>= getPath . (<> "/BigWig/Subcluster/")
        seqIndex <- getGenomeIndex
        let output = dir <> B.unpack nm <> ".bw"
        liftIO $ do
            chrSize <- withGenome seqIndex $ return . getChrSizes
            bedToBigWig output chrSize fl
        |] $ return ()
    ["Subcluster_Merge_Tags"] ~> "Subcluster_Make_BigWig"

 --------------------------------------------------------------------------------
-- Make cell by peak matrix
--------------------------------------------------------------------------------
    nodePar "Call_Peaks" [| \x -> do
        opts <- asks _scatacseq_callpeak_opts
        findPeaks "/Feature/Peak/Subcluster/" opts x
        |] $ return ()
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
    node "Cluster_Peak_Mat" [| \(mats, cls) ->
        subMatrix "/Feature/Peak/Cluster/" mats $ cls^.replicates._2.files
        |] $ return ()
    ["Make_Peak_Mat", "Merged_Cluster"] ~> "Cluster_Peak_Mat"

    node "Compute_Cluster_Peak_RAS" [| computePeakRAS "/Feature/Peak/Cluster/" |] $ return ()
    ["Merge_Peaks", "Cluster_Peak_Mat"] ~> "Compute_Cluster_Peak_RAS"
    node "Cluster_Diff_Peak" [| specificPeaks "/Diff/Peak/Cluster/" |] $ return ()
    ["Call_Peaks_Cluster", "Compute_Cluster_Peak_RAS"] ~> "Cluster_Diff_Peak"


    node "Subcluster_Peak_Mat" [| \(mats, cls) ->
        subMatrix "/Feature/Peak/Subcluster/" mats $ cls^.replicates._2.files
        |] $ return ()
    ["Make_Peak_Mat", "Combine_Clusters"] ~> "Subcluster_Peak_Mat"

    node "Compute_Peak_RAS" [| computePeakRAS "/Feature/Peak/Subcluster/" |] $ return ()
    ["Merge_Peaks", "Subcluster_Peak_Mat"] ~> "Compute_Peak_RAS"
    node "Subcluster_Diff_Peak" [| specificPeaks "/Diff/Peak/Subcluster/" |] $ return ()
    ["Call_Peaks", "Compute_Peak_RAS"] ~> "Subcluster_Diff_Peak"

    {-
    node "Subcluster_Diff_Peak_Prep" [| \(pkList, xs, peaks, ref) -> return $ case ref of
        Nothing -> []
        Just ref' -> 
            let input = flip map xs $ \x ->
                    let pk = fromJust $ lookup (B.pack $ T.unpack $ x^.eid) peaks
                    in x & replicates.traverse.files %~ (\f -> (f, pk))
            in zip3 (repeat $ fromJust pkList) input $ repeat ref'
        |] $ return ()
    ["Merge_Peaks", "Subcluster_Peak_Mat", "Call_Peaks", "Make_Ref_Peak_Mat"] ~> "Subcluster_Diff_Peak_Prep"
    nodePar "Subcluster_Diff_Peak" [| diffPeaks "/Diff/Peak/Subcluster/" |] $ return ()
    path ["Subcluster_Diff_Peak_Prep", "Subcluster_Diff_Peak"]
    -}

--------------------------------------------------------------------------------
-- Make gene matrix
--------------------------------------------------------------------------------
    node "Cluster_Transcript_Mat" [| \(mats, cls) -> 
        subMatrix "/Feature/Gene/Cluster/" mats $ cls^.replicates._2.files
        |] $ return ()
    ["Pre_Make_Transcript_Mat", "Merged_Cluster"] ~> "Cluster_Transcript_Mat"
    nodePar "Compute_Transcript_RAS" [| computeGeneRAS "/Feature/Gene/Cluster/" |] $ return ()
    ["Cluster_Transcript_Mat"] ~> "Compute_Transcript_RAS"
    node "Gene_Acc" [| mkExprTable "/Feature/Gene/Cluster/" |] $ return ()
    ["Pre_Get_Promoters", "Compute_Transcript_RAS"] ~> "Gene_Acc"

    node "Subcluster_Transcript_Mat" [| \(mats, cls) -> 
        subMatrix "/Feature/Gene/Subcluster/" mats $ cls^.replicates._2.files
        |] $ return ()
    ["Pre_Make_Transcript_Mat", "Combine_Clusters"] ~> "Subcluster_Transcript_Mat"
    nodePar "Compute_Subcluster_Transcript_RAS" [| computeGeneRAS "/Feature/Gene/Subcluster/" |] $ return ()
    ["Subcluster_Transcript_Mat"] ~> "Compute_Subcluster_Transcript_RAS"
    node "Subcluster_Gene_Acc" [| mkExprTable "/Feature/Gene/Subcluster/" |] $ return ()
    ["Pre_Get_Promoters", "Compute_Subcluster_Transcript_RAS"] ~> "Subcluster_Gene_Acc"


--------------------------------------------------------------------------------
-- Call CRE interactions
--------------------------------------------------------------------------------
    -- Motif finding
    node "Find_TFBS_Prep" [| findMotifsPre 1e-5 |] $ return ()
    nodePar "Find_TFBS" 'findMotifs $ return ()
    path ["Merge_Peaks", "Find_TFBS_Prep", "Find_TFBS"]
