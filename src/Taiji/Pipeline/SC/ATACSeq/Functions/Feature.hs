{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE GADTs #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Feature
    ( module Taiji.Pipeline.SC.ATACSeq.Functions.Feature.Window
    , module Taiji.Pipeline.SC.ATACSeq.Functions.Feature.Peak
    , module Taiji.Pipeline.SC.ATACSeq.Functions.Feature.Gene
    , module Taiji.Pipeline.SC.ATACSeq.Functions.Feature.Motif
    , dropFeatures
    ) where

import Data.Conduit.Zlib (multiple, ungzip)
import qualified Data.Vector.Unboxed as U
import qualified Data.ByteString.Char8 as B
import Bio.Utils.Functions (scale)
import Control.DeepSeq (force)
import Data.Binary (encodeFile)

import Taiji.Utils.Matrix
import Taiji.Prelude
import Taiji.Pipeline.SC.ATACSeq.Types
import Taiji.Pipeline.SC.ATACSeq.Functions.Feature.Window
import Taiji.Pipeline.SC.ATACSeq.Functions.Feature.Peak
import Taiji.Pipeline.SC.ATACSeq.Functions.Feature.Gene
import Taiji.Pipeline.SC.ATACSeq.Functions.Feature.Motif

dropFeatures :: SCATACSeqConfig config
             => [ SCATACSeq S ( File '[RowName, Gzip] 'Tsv
                              , File '[ColumnName, Gzip] 'Tsv
                              , File '[Gzip] 'Matrix) ]
             -> ReaderT config IO (Maybe FilePath)
dropFeatures input = do
    dir <- asks ((<> "/Spectral/") . _scatacseq_output_dir) >>= getPath
    let output = dir <> "dropped_feats.bin"
    liftIO $ do
        counts <- runConduit $ yieldMany input .| mapMC readCounts .| foldlC f Nothing
        case counts of
            Nothing -> return Nothing
            Just c -> do
                let (zeros, nonzeros) = U.partition ((==0) . snd) $
                        U.zip (U.enumFromN 0 (U.length c)) c
                    (i, v) = U.unzip nonzeros
                    idx = U.toList $ fst $ U.unzip $ U.filter ((>1.65) . snd) $
                        U.zip i $ scale $ U.map fromIntegral v :: [Int]
                encodeFile output $ idx <> U.toList (fst $ U.unzip zeros)
                return $ Just output
  where
    readCounts x = do
        hasCounts <- fmap (fromMaybe False . fmap ((==2) . length . B.split '\t')) $
            runResourceT $ runConduit $ sourceFile columnName .| multiple ungzip .|
                linesUnboundedAsciiC .| headC
        if hasCounts
            then runResourceT $ runConduit $ sourceFile columnName .| multiple ungzip .|
                linesUnboundedAsciiC .| mapC (readInt . last . B.split '\t') .| sinkVector
            else mkSpMatrix readInt mat >>= colSum
      where
        columnName = x^.replicates._2.files._2.location
        mat = x^.replicates._2.files._3.location 
    f Nothing x = force $ Just x
    f (Just x') x = force $ Just $ U.zipWith (+) x' x