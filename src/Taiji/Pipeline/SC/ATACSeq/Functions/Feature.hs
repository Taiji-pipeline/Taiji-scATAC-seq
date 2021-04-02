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
    ) where

import Taiji.Pipeline.SC.ATACSeq.Functions.Feature.Window
import Taiji.Pipeline.SC.ATACSeq.Functions.Feature.Peak
import Taiji.Pipeline.SC.ATACSeq.Functions.Feature.Gene
import Taiji.Pipeline.SC.ATACSeq.Functions.Feature.Motif

{-
readMatrix :: SCATACSeq S (File tags 'Other) -> IO (SpMatrix Int)
readMatrix input = do
    mat <- mkSpMatrix readInt $ input^.replicates._2.files.location
    let prefix = B.pack $ T.unpack (mat^.eid) <> "_" <> show (mat^.replicates._1) <> "+"
    return $ mapRows (first (prefix <>)) mat
{-# INLINE readMatrix #-}

selectFeatures :: [SCATACSeq S (File tags 'Other)] -> 
selectFeatures =
    -}