{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
{-# LANGUAGE LambdaCase #-}
module Taiji.Pipeline.SC.ATACSeq (builder) where

import           Control.Workflow
import qualified Data.Text as T
import qualified Data.ByteString.Char8 as B
import           Bio.Seq.IO (withGenome, getChrSizes)
import           Bio.Data.Bed (readBed, streamBedGzip, sinkFileBedGzip, BED3)
import Data.List.Ordered (nubSort)
import Data.Binary

import           Taiji.Prelude
import           Taiji.Utils
import           Taiji.Pipeline.SC.ATACSeq.Functions
import Taiji.Pipeline.SC.ATACSeq.Types

getFeatures :: SCATACSeqConfig config
            => ( [SCATACSeq S [(B.ByteString, File '[Gzip] 'NarrowPeak)]]
               , [SCATACSeq S (File tags 'Bed, File tags 'Bed, Int)] )
            -> ReaderT config IO (Maybe (File '[Gzip] 'NarrowPeak))
getFeatures ([], []) = return Nothing
getFeatures (peaks, windows) = asks _scatacseq_cluster_by_window >>= \case
    True -> do
        dir <- asks _scatacseq_output_dir >>= getPath . (<> asDir tmp)
        let output = dir <> "/merged.window.bed.gz" 
            f acc x = do
                regions <- liftIO $ runResourceT $ runConduit $
                    streamBedGzip (x^.replicates._2.files._2.location) .|
                    sinkList
                return $ nubSort $ acc ++ (regions :: [BED3])
        beds <- foldM f [] windows
        liftIO $ runResourceT $ runConduit $ yieldMany beds .| sinkFileBedGzip output
        return $ Just $ location .~ output $ emptyFile
    _ -> mergePeaks tmp $ flip concatMap peaks $ \pk -> pk^.replicates._2.files
  where
    tmp = "/temp/Pre/Peak/"

-- | The basic analysis.
basicAnalysis :: Builder ()
basicAnalysis = do
    node "Read_Input" 'readInput $
        doc .= "Read input data information."
    nodePar "Download_Data" 'download $ doc .= "Download data."
    node "Make_Index" 'mkIndices $ doc .= "Generate the genome index."

    uNode "Get_Fastq" [| return . getFastq |]
    nodePar "Demultiplex" 'demultiplex $ return ()
    path ["Read_Input", "Download_Data", "Make_Index", "Get_Fastq", "Demultiplex"]

    uNode "Get_Fastq_Demulti"  [| \(input, x) -> return $ getDemultiFastq input <> x |]
    nodePar "Align" 'tagAlign $ do
        nCore .= 8
        doc .= "Read alignment using BWA. The default parameters are: " <>
            "bwa mem -M -k 32."
    ["Download_Data", "Demultiplex"] ~> "Get_Fastq_Demulti"
    path ["Get_Fastq_Demulti", "Align"]

    uNode "Filter_Bam_Prep" [| \(input, x) -> return $ getBamUnsorted input ++ x |]
    nodePar "Filter_Bam" 'filterNameSortBam $ do
        nCore .= 2
        doc .= "Remove low quality tags using: samtools -F 0x70c -q 30"
    ["Make_Index", "Align"] ~> "Filter_Bam_Prep"
    ["Filter_Bam_Prep"] ~> "Filter_Bam"

    uNode "Get_Bam" [| \(input, x) -> return $ getBam input ++ x |]
    ["Make_Index", "Filter_Bam"] ~> "Get_Bam"

    nodePar "Remove_Duplicates" 'deDuplicates $ return ()
    path ["Get_Bam", "Remove_Duplicates"]

    uNode "Get_Bed" [| \(input, x) -> return $ getSortedBed input ++ x |]
    ["Make_Index", "Remove_Duplicates"] ~> "Get_Bed"

    nodePar "Run_QC" 'getQCMetric $ return ()
    nodePar "Filter_Cell" 'filterCell $ return ()
    path ["Get_Bed", "Run_QC", "Filter_Cell"]

-- PreClustering and doublet detection
preClustering :: Builder ()
preClustering = do
    path ["Filter_Cell", "Pre_Get_Windows"]
    ["Filter_Cell", "Pre_Cluster"] ~> "Pre_Extract_Tags_Prep"
    ["Pre_Make_Peak_Mat", "Run_QC"] ~> "Pre_Detect_Doublet_Prep"
    ["Filter_Cell", "Pre_Detect_Doublet"] ~> "Pre_Remove_Doublets_Prep"
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
        uNode "Extract_Tags_Prep"  [| return . uncurry zipExp |]
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
        nodePar "Call_Peaks" [| \input -> do
            let dir = printf "/temp/Pre/Peak/%s_rep%d/" (T.unpack $ input^.eid) 
                    (input^.replicates._1)
            input & replicates.traverse.files %%~ mapM (findPeaks dir
                 ( def & mode .~ NoModel (-100) 200
                       & cutoff .~ QValue 0.05
                 ))
            |] $ return ()
        path ["Extract_Tags", "Call_Peaks"]

        uNode "Make_Peak_Mat_Prep" [| \(x, y) -> return $ flip map (zipExp x y) $
            \input -> input & replicates._2.files %~
                ( \((a,_,c), pk) ->(a,pk,c) ) |]
        ["Get_Windows", "Call_Peaks"] ~> "Make_Peak_Mat_Prep"
        nodePar "Make_Peak_Mat" [| \input -> do
            let dir = printf "/temp/Pre/Peak/%s_rep%d/" (T.unpack $ input^.eid) 
                    (input^.replicates._1)
            input' <- input & replicates.traverse.files._2 %%~ fmap fromJust .
                mergePeaks dir
            mkPeakMat "/temp/Pre/Peak/" input'
            |] $ return ()
        path ["Make_Peak_Mat_Prep", "Make_Peak_Mat"]

        -- Doublet detection
        uNode "Detect_Doublet_Prep" [| return . uncurry zipExp |]
        nodePar "Detect_Doublet" 'detectDoublet $ return ()
        path ["Detect_Doublet_Prep", "Detect_Doublet"]
        uNode "Remove_Doublets_Prep" [| return . uncurry zipExp |]
        nodePar "Remove_Doublets" 'removeDoublet $ return ()
        path ["Remove_Doublets_Prep", "Remove_Doublets"]

        -- Make Gene matrix 
        node "Get_Promoters" [| \input -> if null input
            then return Nothing
            else Just <$> writePromoters PromoterOnly
            |] $ doc .= "Get the list of promoters from the annotation file."
        ["Remove_Doublets"] ~> "Get_Promoters"
        uNode "Make_Gene_Mat_Prep" [| \(xs, genes) -> return $ zip xs $ repeat $ fromJust genes |]
        ["Remove_Doublets", "Get_Promoters"] ~> "Make_Gene_Mat_Prep"
        nodePar "Make_Gene_Mat" [| mkCellByGene "/temp/Pre/Gene/" |] $
            doc .= "Create cell by transcript matrix for each sample."
        node "Merge_Gene_Mat" [| \mats -> if null mats
            then return Nothing
            else do
                dir <- asks _scatacseq_output_dir >>= getPath . (<> (asDir "/temp/Pre/Gene/"))
                let output = dir <> "Merged_cell_by_gene.mat.gz"
                liftIO $ concatMatrix output $ flip map mats $ \mat ->
                    ( Just $ B.pack $ T.unpack (mat^.eid) <> "_" <> show (mat^.replicates._1)
                    , mat^.replicates._2.files.location )
                return $ Just $ (head mats & eid .~ "Merged") &
                    replicates._2.files.location .~ output
            |] $ return ()
        path ["Make_Gene_Mat_Prep", "Make_Gene_Mat", "Merge_Gene_Mat"]
 
        -- Make Peak matrix
        node "Get_Peak_List" 'getFeatures $ return ()
        ["Call_Peaks", "Get_Windows"] ~> "Get_Peak_List"
        uNode "Make_Feat_Mat_Prep" [| \(bed, pk) -> return $
            flip map (zip bed $ repeat $ fromJust pk) $ \(x, p) ->
                x & replicates.traverse.files %~ (\(a,b) -> (a,p,b))
            |]
        nodePar "Make_Feat_Mat" [| mkPeakMat "/temp/Pre/Feature/" |] $ return ()
        ["Remove_Doublets", "Get_Peak_List"] ~> "Make_Feat_Mat_Prep"
        node "Merge_Feat_Mat" [| \mats -> if null mats
            then return Nothing
            else do
                dir <- asks _scatacseq_output_dir >>= getPath . (<> (asDir "/temp/Pre/Feature/"))
                let output = dir <> "Merged_cell_by_peak.mat.gz"
                liftIO $ concatMatrix output $ flip map mats $ \mat ->
                    ( Just $ B.pack $ T.unpack (mat^.eid) <> "_" <> show (mat^.replicates._1)
                    , mat^.replicates._2.files.location )
                return $ Just $ (head mats & eid .~ "Merged") &
                    replicates._2.files.location .~ output
            |] $ return ()
        path ["Make_Feat_Mat_Prep", "Make_Feat_Mat", "Merge_Feat_Mat"]

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
-- Construct KNN
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
    node "Merged_Batch_Correction" [| \case
        Nothing -> return Nothing
        Just x -> Just <$> batchCorrection "/Cluster/" x
        |] $ return ()
    node "Merged_Make_KNN" [| \case
        Nothing -> return Nothing
        Just x -> Just <$> mkKNNGraph "/Cluster/" x
        |] $ nCore .= 4
    path ["Pre_Merge_Feat_Mat", "Merged_Filter_Mat", "Merged_Reduce_Dims",
        "Merged_Batch_Correction", "Merged_Make_KNN"]

--------------------------------------------------------------------------------
-- Selecting parameter
--------------------------------------------------------------------------------
    uNode "Merged_Param_Search_Prep" [| \case
        (Just spectral, Just knn) -> do
            res <- asks _scatacseq_cluster_resolution_list
            optimizer <- asks _scatacseq_cluster_optimizer 
            return $ flip map res $ \r ->
                ( optimizer, r
                , spectral^.replicates._2.files._2.location
                , knn^.replicates._2.files._2.location )
        _ -> return []
        |]
    ["Merged_Reduce_Dims", "Merged_Make_KNN"] ~> "Merged_Param_Search_Prep"
    nodePar "Merged_Param_Search" [| \(optimizer, r, spectral, knn) -> do
        res <- liftIO $ evalClusters optimizer r spectral knn
        return (r, res)
        |] $ nCore .= 5
    path ["Merged_Param_Search_Prep", "Merged_Param_Search"]

    node "Merged_Get_Param" [| \(knn, res) -> case knn of
        Nothing -> return Nothing
        Just knn' -> do
            dir <- asks _scatacseq_output_dir >>= getPath . (<> "/Figure/")
            p <- liftIO $ optimalParam (dir <> "Clustering_parameters.html") res
            asks _scatacseq_cluster_resolution >>= \case
                Nothing -> return $ Just (p, knn')
                Just p' -> return $ Just (p', knn')
        |] $ return ()
    ["Merged_Make_KNN", "Merged_Param_Search"] ~> "Merged_Get_Param"
    node "Merged_Cluster" [| \case
        Nothing -> return Nothing
        Just (res, input) -> do
            optimizer <- asks _scatacseq_cluster_optimizer 
            Just <$> clustering "/Cluster/" res optimizer input
        |] $ return ()
    path ["Merged_Get_Param", "Merged_Cluster"]

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
    uNode "Extract_Tags_Prep" [| \(x,y) -> return $ zip x $ repeat $ fromJust y |]
    extractTags "/Bed/Cluster/"
    ["Pre_Remove_Doublets", "Merged_Cluster"] ~> "Extract_Tags_Prep"
    ["Extract_Tags_Prep"] ~> "Extract_Tags"

    nodePar "Make_BigWig" [| \(nm, fl) -> do
        dir <- asks _scatacseq_output_dir >>= getPath . (<> "/BigWig/Cluster/")
        seqIndex <- asks ( fromMaybe (error "Genome index file was not specified!") .
            _scatacseq_genome_index )
        tmpdir <- fromMaybe "./" <$> asks _scatacseq_tmp_dir
        let output = dir <> B.unpack nm <> ".bw"
        blackRegions <- asks _scatacseq_blacklist >>= \case
            Nothing -> return []
            Just blacklist -> liftIO $ readBed blacklist
        liftIO $ do
            chrSize <- withGenome seqIndex $ return . getChrSizes
            bedToBigWig output chrSize blackRegions tmpdir fl
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

    uNode "Make_Peak_Mat_Prep" [| \(bed, pk) -> return $
        flip map (zip bed $ repeat $ fromJust pk) $ \(x, p) ->
            x & replicates.traverse.files %~ (\(a,b) -> (a,p,b))
        |]
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
-- Subclustering 
--------------------------------------------------------------------------------
    uNode "Subcluster_Reduce_Dims_Prep" [| \xs -> 
        let f x = do 
                nCells <- liftIO $ _num_row <$>
                    mkSpMatrix id (x^.replicates._2.files.location)
                return $ nCells >= 1000
        in filterM f xs
        |]
    nodePar "Subcluster_Reduce_Dims" 'subSpectral $ return ()
    nodePar "Subcluster_Make_KNN" [| mkKNNGraph "/Subcluster/KNN/" |] $ return ()
    path ["Cluster_Peak_Mat", "Subcluster_Reduce_Dims_Prep", "Subcluster_Reduce_Dims", "Subcluster_Make_KNN"]

    uNode "Subcluster_Param_Search_Prep" [| \(sp, knn) -> do
        res <- asks _scatacseq_cluster_resolution_list
        optimizer <- asks _scatacseq_cluster_optimizer 
        return $ flip map (zipExp sp knn) $ \e -> e & replicates.traversed.files %~
            ( \(x,y) -> flip map res $ \r -> (optimizer, r, x^._2.location, y^._2.location) )
        |]
    nodePar "Subcluster_Param_Search" [| \input -> input & replicates.traversed.files %%~
        liftIO . ( \xs -> do
            res <- forM xs $ \(optimizer, r, spectral, knn) -> do
                res <- liftIO $ evalClusters optimizer r spectral knn
                return (r, res)
            return (head xs ^. _4, res)
            )
        |] $ nCore .= 5
    nodePar "Subcluster_Get_Param" [| \input -> do
        dir <- asks _scatacseq_output_dir >>= getPath . (<> "/Subcluster/")
        let output = printf "%s/%s_parameters.html" dir (T.unpack $ input^.eid)
        input & replicates.traversed.files %%~ ( \(knn, res) -> do
            p <- liftIO $ optimalParam output res
            return (p, knn)
            )
        |] $ return ()
    ["Subcluster_Reduce_Dims", "Subcluster_Make_KNN"] ~> "Subcluster_Param_Search_Prep"
    path ["Subcluster_Param_Search_Prep", "Subcluster_Param_Search", "Subcluster_Get_Param"]
 


--------------------------------------------------------------------------------
-- Gene accessibility
--------------------------------------------------------------------------------
    node "Gene_Acc" [| mkExprTable "/Feature/Gene/Cluster/" |] $ return ()
    ["Pre_Get_Promoters", "Peak_Acc"] ~> "Gene_Acc"

--------------------------------------------------------------------------------
-- Run ChromVar 
--------------------------------------------------------------------------------
    -- Motif finding
    node "Find_TFBS_Prep" [| findMotifsPre 1e-5 |] $ return ()
    nodePar "Find_TFBS" 'findMotifs $ return ()
    path ["Merge_Peaks", "Find_TFBS_Prep", "Find_TFBS"]

    -- ChromVar
    {-
    node "Make_Motif_Peak_Mat" 'mkMotifMat $ return ()
    ["Merge_Peaks", "Find_TFBS"] ~> "Make_Motif_Peak_Mat"

    node "ChromVar_Pre" 'preChromVar $ return ()
    ["Make_Motif_Peak_Mat", "Merge_Peaks", "Cluster_Peak_Mat"] ~> "ChromVar_Pre"
    nodePar "ChromVar" 'runChromVar $ return ()
    ["ChromVar_Pre"] ~> "ChromVar"
    -}