{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE GADTs #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Feature
    ( genPeakMat
    , mergeFeatMatrix
    , module Taiji.Pipeline.SC.ATACSeq.Functions.Feature.Window
    , module Taiji.Pipeline.SC.ATACSeq.Functions.Feature.Peak
    , module Taiji.Pipeline.SC.ATACSeq.Functions.Feature.Gene
    , module Taiji.Pipeline.SC.ATACSeq.Functions.Feature.Motif
    ) where

import qualified Data.ByteString.Char8 as B
import Data.Conduit.List (groupBy)
import           Control.Workflow
import Bio.Data.Bed.Types
import Bio.Data.Bed
import Data.Singletons.Prelude (Elem)
import           Bio.Pipeline.Utils
import qualified Data.Text as T
import Bio.Seq.IO (withGenome, getChrSizes)

import Taiji.Prelude hiding (groupBy)
import Taiji.Pipeline.SC.ATACSeq.Types
import Taiji.Pipeline.SC.ATACSeq.Functions.Utils
import Taiji.Pipeline.SC.ATACSeq.Functions.Feature.Window
import Taiji.Pipeline.SC.ATACSeq.Functions.Feature.Peak
import Taiji.Pipeline.SC.ATACSeq.Functions.Feature.Gene
import Taiji.Pipeline.SC.ATACSeq.Functions.Feature.Motif

genPeakMat :: FilePath   -- ^ directory
           -> Maybe T.Text     -- ^ namespace
           -> T.Text     -- ^ Input 1
           -> T.Text     -- ^ Input 2
           -> Builder ()
genPeakMat prefix nm input1 input2 = case nm of
    Just n -> do
        namespace n builder
        path [input1, n <> "_Call_Peaks"]
        [input2, n <> "_Merge_Peaks"] ~> (n <> "_Make_Peak_Matrix_Prep")
    Nothing -> do
        builder
        path [input1, "Call_Peaks"]
        [input2, "Merge_Peaks"] ~> "Make_Peak_Matrix_Prep"
  where
    builder = do
        nodePar "Call_Peaks" [| findPeaks $ prefix <> "/Peaks/" |] $ return ()
        node "Merge_Peaks" [| mergePeaks prefix |] $ return ()
        path ["Call_Peaks", "Merge_Peaks"]

        node "Make_Peak_Matrix_Prep" [| \(exps, pk) -> return $ flip map exps $
            \e -> e & replicates._2.files %~ (\(a,_,c) -> (a,fromJust pk,c))
            |] $ return ()
        nodePar "Make_Peak_Matrix" [| mkPeakMat prefix |] $ return ()

        path ["Make_Peak_Matrix_Prep", "Make_Peak_Matrix"]
        node "Merge_Peak_Matrix_Prep" [| \(pk, exps) -> return $ flip map exps $
            \e -> e & replicates._2.files %~ (\x -> (fromJust pk,x))
            |]$ return ()
        node "Merge_Peak_Matrix" [| mergeFeatMatrix $ prefix ++ "/merged_cell_by_peak" |] $ return ()
        ["Merge_Peaks", "Make_Peak_Matrix"] ~> "Merge_Peak_Matrix_Prep"
        path ["Merge_Peak_Matrix_Prep", "Merge_Peak_Matrix"]


mergeFeatMatrix :: ( Elem 'Gzip tags1 ~ 'True
                   , Elem 'Gzip tags2 ~ 'True
                   , SCATACSeqConfig config )
                => FilePath
                -> [SCATACSeq S (File tags1 file, File tags2 'Other)]
                -> ReaderT config IO [SCATACSeq S (File '[Gzip] 'Bed, File '[Gzip] 'Other)]
mergeFeatMatrix _ [] = return []
mergeFeatMatrix filename inputs = do
    dir <- asks _scatacseq_output_dir >>= getPath
    let output = dir <> "/" <> filename <> ".mat.gz"
        outputIdx = dir <> "/" <> filename <> ".idx.bed.gz"
    liftIO $ runResourceT $ runConduit $ mergeMatrix inputs' outputIdx .| sinkFile output
    return $ return $ (head inputs & eid .~ "Merged") & replicates._2.files .~
        ( location .~ outputIdx $ emptyFile
        , location .~ output $ emptyFile )
  where
    inputs' = map (\x -> (B.pack $ T.unpack $ x^.eid, x^.replicates._2.files)) inputs

--------------------------------------------------------------------------------
-- Cut site map
--------------------------------------------------------------------------------

-- | Create the cut site map for every cell.
mkCutSiteIndex :: SCATACSeqConfig config
               => SCATACSeq S (File '[NameSorted, Gzip] 'Bed)
               -> ReaderT config IO (SCATACSeq S (File '[] 'Other))
mkCutSiteIndex input = do
    dir <- asks ((<> "/CutSiteIndex") . _scatacseq_output_dir) >>= getPath
    genome <- fromJust <$> asks _scatacseq_genome_index 
    let output = printf "%s/%s_rep%d.csidx" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . (\fl -> do
        chrs <- withGenome genome $ return . map fst . getChrSizes
        mkCutSiteIndex_ output fl chrs
        return $ emptyFile & location .~ output )
{-# INLINE mkCutSiteIndex #-}

mkCutSiteIndex_ :: FilePath
                -> File '[NameSorted, Gzip] 'Bed
                -> [B.ByteString]
                -> IO ()
mkCutSiteIndex_ output input chrs = createCutSiteIndex output chrs $
    streamBedGzip (input^.location) .| groupBy ((==) `on` (^.name)) .|
    mapC (\x -> (fromJust $ head x ^. name, x))
{-# INLINE mkCutSiteIndex_ #-}