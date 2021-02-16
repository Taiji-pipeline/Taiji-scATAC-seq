{-# OPTIONS_GHC -fno-full-laziness #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE LambdaCase #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Internal (demulti) where

import Bio.Data.Fastq
import Data.Conduit.Internal (zipSources, zipSinks)
import qualified Data.ByteString.Char8                as B
import Data.Maybe
import qualified Data.IntMap.Strict as I
import Conduit
import Bio.Pipeline

demulti :: FilePath -> FilePath -> I.IntMap Int -> Int -> FilePath -> FilePath -> FilePath -> IO ()
demulti out1 out2 bcMap k fqidx fq1 fq2 = do
    _ <- runResourceT $ runConduit $ zipSources
        (streamFastqGzip fqidx .| mapC (dnaToInt . B.take k . fastqSeq))
        (zipSources (streamFastqGzip fq1) (streamFastqGzip fq2)) .|
        f bcMap .| zipSinks (mapC fst .| sinkFastqGzip out1) (mapC snd .| sinkFastqGzip out2)
    return ()
  where
    f :: Monad m => I.IntMap Int -> ConduitT (Int, (Fastq, Fastq)) (Fastq, Fastq) m ()
    f bcMap = concatMapC $ \(bc, (fq1, fq2)) -> case I.lookup bc bcMap of
        Nothing -> Nothing
        Just bc' -> Just (addBarcode bc' fq1, addBarcode bc' fq2)
      where
        addBarcode bc fq = fq{fastqSeqId = B.concat [intToDna bc, ":", fastqSeqId fq]}