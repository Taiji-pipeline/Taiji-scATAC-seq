{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE GADTs #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Feature
    ( mergeFeatMatrix
    , module Taiji.Pipeline.SC.ATACSeq.Functions.Feature.Window
    , module Taiji.Pipeline.SC.ATACSeq.Functions.Feature.Peak
    , module Taiji.Pipeline.SC.ATACSeq.Functions.Feature.Gene
    , module Taiji.Pipeline.SC.ATACSeq.Functions.Feature.Motif
    ) where

import qualified Data.ByteString.Char8 as B
import Data.Conduit.List (groupBy)
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

mergeFeatMatrix :: ( Elem 'Gzip tags1 ~ 'True
                   , Elem 'Gzip tags2 ~ 'True
                   , SCATACSeqConfig config )
                => FilePath
                -> [SCATACSeq S (File tags1 file, File tags2 'Other)]
                -> ReaderT config IO [SCATACSeq S (File '[Gzip] 'Other)]
mergeFeatMatrix _ [] = return []
mergeFeatMatrix filename inputs = do
    dir <- asks ((<> "/Feature") . _scatacseq_output_dir) >>= getPath
    let output = dir ++ "/" ++ filename
    liftIO $ runResourceT $ runConduit $ mergeMatrix inputs' .| sinkFile output
    return $ return $ (head inputs & eid .~ "Merged") & replicates._2.files .~
        (location .~ output $ emptyFile)
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