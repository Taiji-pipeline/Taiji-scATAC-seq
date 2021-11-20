{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
{-# LANGUAGE LambdaCase #-}
module Taiji.Pipeline.SC.ATACSeq (builder) where

import           Control.Workflow
import qualified Data.Text as T
import qualified Data.ByteString.Char8 as B
import qualified Data.Map.Strict as Map
import           Bio.Seq.IO (withGenome, getChrSizes)
import           Bio.Data.Bed (readBed, streamBedGzip, sinkFileBedGzip, BED3, mergeSortedBed, splitBedBySizeLeft, mergeBed)
import Data.List.Ordered (nubSort)
import Data.Binary
import Control.DeepSeq (($!!))

import           Taiji.Prelude
import           Taiji.Utils
import           Taiji.Pipeline.SC.ATACSeq.Functions
import Taiji.Pipeline.SC.ATACSeq.Types

getFeatures :: SCATACSeqConfig config
            => ( [SCATACSeq S [(B.ByteString, File '[Gzip] 'NarrowPeak)]]
               , [SCATACSeq S (a, File tags 'Bed, Int)] )
            -> ReaderT config IO (Maybe (File '[Gzip] 'NarrowPeak))
getFeatures ([], []) = return Nothing
getFeatures (peaks, windows) = do
    dir <- asks _scatacseq_output_dir >>= getPath . (<> asDir tmp)
    ws <- fromIntegral <$> asks _scatacseq_window_size
    asks _scatacseq_cluster_by_window >>= \case
        True -> do
            let output = dir <> "/merged.window.bed.gz" 
                f acc x = do
                    regions <- liftIO $ runResourceT $ runConduit $
                        streamBedGzip (x^.replicates._2.files._2.location) .|
                        mergeSortedBed .| sinkList
                    return $!! runIdentity $ runConduit $
                        mergeBed (acc <> regions :: [BED3]) .| sinkList
            beds <- foldM f [] windows
            liftIO $ runResourceT $ runConduit $ yieldMany beds .|
                concatMapC (splitBedBySizeLeft ws) .|
                sinkFileBedGzip output
            return $ Just $ location .~ output $ emptyFile
        _ -> mergePeaks dir $ concatMap (^.replicates._2.files) peaks
  where
    tmp = "/temp/Peak/"

-- | The basic analysis.
basicAnalysis :: Builder ()
basicAnalysis = do
    node "Read_Input" 'readInput $
        doc .= "Read input data information."
    nodePar "Download_Data" 'download $ do
        doc .= "Download data."
        nCore .= 2
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

    uNode "Get_Bed" [| return . getBedFiles |]
    nodePar "Sort_Bed" 'sortBedFile $ return ()
    path ["Make_Index", "Get_Bed", "Sort_Bed"]

    uNode "Get_Sorted_Bed" [| \(x, y) -> return $ x ++ y |]
    ["Sort_Bed", "Remove_Duplicates"] ~> "Get_Sorted_Bed"

    nodePar "Run_QC" 'getQCMetric $ nCore .= 2
    path ["Get_Sorted_Bed", "Run_QC"]

-- PreClustering and doublet detection
preClustering :: Builder ()
preClustering = do
    path ["Run_QC", "Pre_Get_Windows"]
    ["Run_QC", "Pre_Cluster"] ~> "Pre_Call_Peaks_Prep"
    ["Run_QC", "Pre_Call_Peaks"] ~> "Pre_Detect_Doublet_Prep"
    namespace "Pre" $ do
        nodePar "Get_Windows" [| getWindows "/Feature/Window/" |] $ return ()

        -- Clustering in each sample
        nodePar "Cluster" [| \input -> do
            dir <- asks _scatacseq_tmp_dir
            let prefix = "/temp/Cluster/"
            withRunInIO $ \runInIO -> withTempDir dir $ \tmp -> runInIO $ do
                mkWindowMat tmp input >>= filterMatrix tmp >>=
                    spectral tmp Nothing >>= mkKNNGraph tmp >>=
                    clustering prefix 1 RBConfigurationWeighted
            |] $ return ()
        path ["Get_Windows", "Cluster"]

        uNode "Call_Peaks_Prep"  [| return . uncurry zipExp |]
        nodePar "Call_Peaks" [| \input -> do
            dir <- asks _scatacseq_tmp_dir
            passedQC <- getQCFunction
            withRunInIO $ \runInIO -> withTempDir dir $ \tmp -> runInIO $ do
                let outDir = printf "/temp/Peak/%s_rep%d/" (T.unpack $ input^.eid) 
                        (input^.replicates._1)
                input & replicates.traverse.files %%~ ( \(tags, cl) -> do
                    clusters <- liftIO $ decodeFile $ cl^.location
                    bedFls <- liftIO $ runResourceT $ runConduit $ streamTN5Insertion passedQC tags .|
                        splitBedByCluster' tmp clusters
                    forM bedFls $ findPeaks outDir
                        ( def & mode .~ NoModel (-100) 200
                            & cutoff .~ QValue 0.05
                            & tmpDir .~ tmp
                        )
                    )
            |] $ return ()
        path ["Call_Peaks_Prep", "Call_Peaks"]

        -- Doublet detection
        uNode "Detect_Doublet_Prep" [| \(x, y) -> return $ zipExp x y |]
        nodePar "Detect_Doublet" [| \input -> do
            dir <- asks _scatacseq_tmp_dir
            withRunInIO $ \runInIO -> withTempDir dir $ \tmp -> runInIO $ do
                input' <- input & replicates.traverse.files._2 %%~ fmap fromJust .
                    mergePeaks tmp
                mat <- mkPeakMat tmp input'
                detectDoublet $ mat & replicates.traverse.files %~ ( \x ->
                    (x, input'^.replicates._2.files._1) )
            |] $ return ()
        path ["Detect_Doublet_Prep", "Detect_Doublet"]

        -- Make Peak matrix
        node "Get_Peak_List" 'getFeatures $ return ()
        ["Call_Peaks", "Get_Windows"] ~> "Get_Peak_List"
        uNode "Make_Feat_Mat_Prep" [| \(bed, pk) -> return $
            flip map (zip bed $ repeat $ fromJust pk) $ \(x, p) ->
                x & replicates.traverse.files %~ (\a -> (a,p))
            |]
        nodePar "Make_Feat_Mat" [| \input -> do
            dir <- asks ((<> "/Feature/Sample/") . _scatacseq_output_dir) >>= getPath
            mkFeatMat dir input
            |] $ return ()
        ["Detect_Doublet", "Get_Peak_List"] ~> "Make_Feat_Mat_Prep"
        path ["Make_Feat_Mat_Prep", "Make_Feat_Mat"]

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
-- Gene matrix
--------------------------------------------------------------------------------
    node "Get_Promoters" [| \input -> if null input
        then return Nothing
        else Just <$> writePromoters PromoterOnly
        |] $ doc .= "Get the list of promoters from the annotation file."
    ["Pre_Detect_Doublet"] ~> "Get_Promoters"
    uNode "Make_Gene_Mat_Prep" [| \(xs, genes) -> return $ zip xs $ repeat $ fromJust genes |]
    nodePar "Make_Gene_Mat" [| mkCellByGene "/Feature/Gene/Sample/" |] $
        doc .= "Create cell by transcript matrix for each sample."
    ["Pre_Detect_Doublet", "Get_Promoters"] ~> "Make_Gene_Mat_Prep"
    path ["Make_Gene_Mat_Prep", "Make_Gene_Mat"]

--------------------------------------------------------------------------------
-- Construct KNN
--------------------------------------------------------------------------------
    uNode "Get_Feat_Mat" [| \(input, mats) -> case getMatrix input of
        [] -> return mats
        x -> return x
        |]
    ["Read_Input", "Pre_Make_Feat_Mat"] ~> "Get_Feat_Mat"
    node "Merged_Feature_Selection" 'dropFeatures $ return ()
    path ["Get_Feat_Mat", "Merged_Feature_Selection"]

    node "Merged_Spectral" [| getSpectral 35000 |] $ nCore .= 4
    ["Get_Feat_Mat", "Merged_Feature_Selection"] ~> "Merged_Spectral"

    uNode "Merged_Nystrom_Prep" [| \(input, model, feat) -> case model of
        Just (Right m) -> liftIO $ do
            input' <- chunksInput $ map (\x -> x^.replicates._2.files._3) input
            return $ zip3 input' (repeat m) $ repeat $ fromJust feat
        _ -> return []
        |]
    nodePar "Merged_Nystrom" 'nystromExtend $ return ()
    ["Get_Feat_Mat", "Merged_Spectral", "Merged_Feature_Selection"] ~> "Merged_Nystrom_Prep"
    path ["Merged_Nystrom_Prep", "Merged_Nystrom"]

    node "Merged_Reduce_Dims" 'mergeResults $ return ()
    ["Merged_Spectral", "Get_Feat_Mat", "Merged_Nystrom"] ~> "Merged_Reduce_Dims"

    node "Merged_Batch_Correction" [| \case
        Nothing -> return Nothing
        Just x -> Just <$> batchCorrection "/Cluster/" x
        |] $ return ()
    node "Merged_Make_KNN" [| \case
        Nothing -> return Nothing
        Just x -> do
            dir <- asks ((<> "/Cluster/") . _scatacseq_output_dir) >>= getPath
            Just <$> mkKNNGraph dir x
        |] $ nCore .= 4
    path ["Merged_Reduce_Dims", "Merged_Batch_Correction", "Merged_Make_KNN"]

--------------------------------------------------------------------------------
-- Clustering
--------------------------------------------------------------------------------
    uNode "Merged_Param_Search_Prep" [| \case
        Just knn -> do
            ress <- asks _scatacseq_cluster_resolution_list
            res <- maybe [] return <$> asks _scatacseq_cluster_resolution
            optimizer <- asks _scatacseq_cluster_optimizer 
            return $ flip map (nubSort $ res <> ress) $ \r ->
                ( optimizer, r
                , knn^.replicates._2.files._2.location )
        _ -> return []
        |]
    nodePar "Merged_Param_Search" [| \(optimizer, r, knn) -> do
        dir <- asks _scatacseq_output_dir >>= getPath . (<> asDir ("/Cluster/Params/" <> show r))
        res <- liftIO $ evalClusters dir optimizer r knn
        return (r, res)
        |] $ return ()
    path ["Merged_Make_KNN", "Merged_Param_Search_Prep", "Merged_Param_Search"]

    uNode "Merged_Cluster_Metric_Prep" [| \(spec, x) -> case spec of
        Nothing -> return []
        Just fl -> return $ flip map x $ \(r, y) -> (r, y, fl^.replicates._2.files._2.location)
        |]
    nodePar "Merged_Cluster_Metric" [| \(res, cl, spec) -> do
        r <- liftIO $ computeClusterMetrics cl spec
        return (res, r)
        |] $ return ()
    node "Merged_Cluster" 'plotClusters $ return ()
    ["Merged_Reduce_Dims", "Merged_Param_Search"] ~> "Merged_Cluster_Metric_Prep"
    path ["Merged_Cluster_Metric_Prep", "Merged_Cluster_Metric"]
    ["Merged_Cluster_Metric", "Merged_Param_Search", "Merged_Make_KNN"] ~> "Merged_Cluster"



--------------------------------------------------------------------------------
-- Subclustering 
--------------------------------------------------------------------------------
    uNode "Subcluster_Get_Features_Prep" [| \case
        (Just clFl, x) -> asks _scatacseq_do_subclustering >>= \case
            False -> return []
            True -> liftIO $ do
                let fl = clFl^.replicates._2.files
                    minimalCells = 100
                clusters <- decodeFile $ fl^.location
                let clNames = map (T.pack . B.unpack . fst) $
                        filter ((>=minimalCells) . snd) $ flip map clusters $ \cl ->
                            (_cluster_name cl, length $ _cluster_member cl)
                return $ flip map clNames $ \nm ->
                    (eid .~ nm $ clFl, x)
        _ -> return []
        |]
    nodePar "Subcluster_Get_Features" 'subsetFeatMat $ return ()
    ["Merged_Cluster", "Get_Feat_Mat"] ~> "Subcluster_Get_Features_Prep"
    path ["Subcluster_Get_Features_Prep", "Subcluster_Get_Features"]
 
    nodePar "Subcluster_Reduce_Dims" 'subSpectral $ return ()
    nodePar "Subcluster_Make_KNN" [| \input -> do
        dir <- asks ((<> "/Subcluster/KNN/") . _scatacseq_output_dir) >>= getPath
        mkKNNGraph dir input
        |] $ return ()
    path ["Subcluster_Get_Features", "Subcluster_Reduce_Dims", "Subcluster_Make_KNN"]
    uNode "Subcluster_Param_Search_Prep" [| \(sp, knn) -> do
        rs <- asks _scatacseq_cluster_resolution_list
        optimizer <- asks _scatacseq_cluster_optimizer 
        fmap (concatMap split) $ forM (zipExp sp knn) $ \e ->
            e & replicates.traversed.files %%~ ( \(x,y) -> do
                r <- Map.lookup (e^.eid) . fromMaybe Map.empty <$>
                    asks _scatacseq_subcluster_resolution 
                let res = sort $ case r of
                        Nothing -> rs
                        Just c -> c : rs
                return $ flip map res $ \r' -> (optimizer, r', x^._2.location, y)
                )
        |]
    nodePar "Subcluster_Param_Search" [| \input -> input & replicates.traversed.files %%~
        ( \(optimizer, r, spec, knn) -> do
            let i = "/Subcluster/Params/" <> T.unpack (input^.eid) <> "/" <> show r <> "/"
            dir <- asks _scatacseq_output_dir >>= getPath . (<> asDir i)
            liftIO $ do
                grs <- evalClusters dir optimizer r (knn^._2.location)
                metric <- computeClusterMetrics grs spec
                return (r, metric, grs, knn)
        )
        |] $ return ()
    uNode "Subcluster_Cluster_Prep" [| \input -> return $
        flip map (concatMap split $ mergeExp input) $ \e -> 
            let (a,b,c,_) = unzip4 $ e^.replicates._2.files
            in (zip a b, zip a c, e & replicates.traversed.files %~ (^._4) . head)
        |]
    nodePar "Subcluster_Cluster" 'plotSubclusters $ return ()
    ["Subcluster_Reduce_Dims", "Subcluster_Make_KNN"] ~> "Subcluster_Param_Search_Prep"
    path ["Subcluster_Param_Search_Prep", "Subcluster_Param_Search",
        "Subcluster_Cluster_Prep", "Subcluster_Cluster"]

    node "Combine_Clusters" 'combineClusters $ return ()
    node "Cluster_Viz" 'vizCluster $ nCore .= 8
    ["Merged_Cluster", "Subcluster_Cluster"] ~> "Combine_Clusters"
    ["Combine_Clusters", "Merged_Batch_Correction", "QC"] ~> "Cluster_Viz"

--------------------------------------------------------------------------------
-- Compute gene-level accessibility
--------------------------------------------------------------------------------
    uNode "Gene_Count_Prep" [| \case
        (mats, Just promoter, Just clFl) -> return $ zip3 mats (repeat promoter) (repeat clFl)
        _ -> return []
        |]
    nodePar "Gene_Count" [| mkExprTable "/temp/Gene/" |] $ return ()
    ["Make_Gene_Mat", "Get_Promoters", "Combine_Clusters"] ~> "Gene_Count_Prep"
    path ["Gene_Count_Prep", "Gene_Count"]
    node "Gene_Acc" [| \input -> do
        dir <- asks _scatacseq_output_dir >>= getPath . (<> "/Feature/Gene/")
        let output = dir <> "gene_acc.tsv"
        combineExprTable output input
        |] $ return ()
    ["Get_Promoters", "Gene_Count"] ~> "Gene_Acc"

    -- Extract tags for each subcluster
    uNode "Extract_Tags_Prep" [| \(x,y) -> return $ zip x $ repeat $ fromJust y |]
    extractTags "/Bed/Cluster/"
    ["Extract_Tags_Prep"] ~> "Extract_Tags"
    ["Pre_Detect_Doublet", "Combine_Clusters"] ~> "Extract_Tags_Prep"

    nodePar "Call_Peaks" [| \x -> do
        opts <- getCallPeakOpt
        findPeaks "/Feature/Peak/Cluster/" opts x
        |] $ return ()
    node "Merge_Peaks" [| \input -> do
        dir <- asks _scatacseq_output_dir >>= getPath . (<> "/Feature/Peak/")
        mergePeaks dir input
        |] $ return ()
    path ["Merge_Tags", "Call_Peaks", "Merge_Peaks"]

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

    uNode "Make_Peak_Mat_Prep" [| \(bed, pk) -> return $
        flip map (zip bed $ repeat $ fromJust pk) $ \(x, p) ->
            x & replicates.traverse.files %~ (\a -> (a,p))
        |]
    nodePar "Make_Peak_Mat" [| \input -> do
        dir <- asks ((<> "/Feature/Peak/Sample/") . _scatacseq_output_dir) >>= getPath
        mkPeakMat dir input
        |] $ return ()
    ["Pre_Detect_Doublet", "Merge_Peaks"] ~> "Make_Peak_Mat_Prep"
    path ["Make_Peak_Mat_Prep", "Make_Peak_Mat"]

    node "Cluster_Peak_Mat" [| \(mats, cl) -> case cl of
        Nothing -> return []
        Just x -> subMatrix "/Feature/Peak/Cluster/" mats $ x^.replicates._2.files
        |] $ return ()
    ["Make_Peak_Mat", "Combine_Clusters"] ~> "Cluster_Peak_Mat"

    node "Cluster_Peak_Acc" [| computePeakRAS "/Feature/Peak/Cluster/" |] $ return ()
    ["Merge_Peaks", "Cluster_Peak_Mat"] ~> "Cluster_Peak_Acc"


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