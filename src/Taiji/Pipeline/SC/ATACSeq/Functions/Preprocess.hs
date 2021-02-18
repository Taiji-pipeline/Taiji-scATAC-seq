{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE GADTs #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Preprocess
    ( readInput
    , download
    , getFastq
    , getDemultiFastq
    , getBamUnsorted
    , getBam
    , getSortedBed
    , demultiplex
    ) where

import           Data.Bifunctor                (bimap)
import Bio.Pipeline
import Data.Either (lefts)
import Bio.Data.Experiment.Parser
import qualified Data.Vector.Unboxed as U
import qualified Data.Text as T
import qualified Data.ByteString.Char8 as B
import Shelly hiding (FilePath)
import Data.Either
import Bio.Data.Fastq (streamFastqGzip)

import Taiji.Prelude
import Taiji.Pipeline.SC.ATACSeq.Types
import Taiji.Pipeline.SC.ATACSeq.Functions.Internal (demulti)

type RAWInput = SCATACSeq N [Either SomeFile (SomeFile, SomeFile)]

readInput :: SCATACSeqConfig config
          => () -> ReaderT config IO [RAWInput]
readInput _ = do
    input <- asks _scatacseq_input
    liftIO $ mkInputReader input "scATAC-seq" SCATACSeq

download :: SCATACSeqConfig config
         => [RAWInput]
         -> ReaderT config IO [RAWInput]
download input = do
    tmp <- fromMaybe "./" <$> asks _scatacseq_tmp_dir
    input & traverse.replicates.traverse.files.traverse %%~ ( \fl -> do
        dir <- asks _scatacseq_output_dir >>= getPath . (<> (asDir "/Download"))
        liftIO $ downloadFiles dir tmp fl )

getFastq :: [RAWInput]
          -> [ SCATACSeq S ( (File '[Gzip] 'Fastq, File '[Gzip] 'Fastq)
                        , File '[Gzip] 'Fastq ) ]
getFastq input = concatMap split $ concatMap split $
    input & mapped.replicates.mapped.files %~ f
  where
    f fls = case partitionEithers fls' of
        ([fqidx], [(fq1, fq2)]) -> [((fromSomeFile fq1, fromSomeFile fq2), fromSomeFile fqidx)]
        _ -> []
      where
        fls' = filter (either (\x -> getFileType x == Fastq) g) fls
        g (x,y) = getFileType x == Fastq && getFileType y == Fastq

getDemultiFastq :: [RAWInput]
         -> [SCATACSeq S (File '[Demultiplexed, Gzip] 'Fastq, File '[Demultiplexed, Gzip] 'Fastq)]
getDemultiFastq input = concatMap split $ concatMap split $
    input & mapped.replicates.mapped.files %~ (\fls -> map (\(x,y) -> (fromSomeFile x, fromSomeFile y)) $
      filter (\(x,y) -> f x && f y) $ rights fls)
  where
    f x = getFileType x == Fastq && x `hasTag` Gzip && x `hasTag` Demultiplexed

getBamUnsorted :: [RAWInput]
               -> [ SCATACSeq S ( Either
               (File '[] 'Bam) (File '[PairedEnd] 'Bam) )]
getBamUnsorted input = concatMap split $ concatMap split $
    input & mapped.replicates.mapped.files %~ mapMaybe f . lefts
  where
    f fl | getFileType fl == Bam && not (fl `hasTag` NameSorted) &&
              fl `hasTag` PairedEnd = Just $ Right $ fromSomeFile fl
         | getFileType fl == Bam && not (fl `hasTag` NameSorted) = Just $ Left $ fromSomeFile fl
         | otherwise = Nothing
 
getBam :: [RAWInput]
       -> [ SCATACSeq S ( Either
          (File '[NameSorted] 'Bam) (File '[NameSorted, PairedEnd] 'Bam) )]
getBam input = concatMap split $ concatMap split $
    input & mapped.replicates.mapped.files %~ mapMaybe f . lefts
  where
    f fl | getFileType fl == Bam && fl `hasTag` NameSorted &&
              fl `hasTag` PairedEnd = Just $ Right $ fromSomeFile fl
         | getFileType fl == Bam && fl `hasTag` NameSorted = Just $ Left $ fromSomeFile fl
         | otherwise = Nothing

getSortedBed :: ( unpair ~ '[NameSorted, Gzip]
                , paired ~ '[NameSorted, PairedEnd, Gzip] )
             => [RAWInput]
             -> [SCATACSeq S (Either (File unpair 'Bed) (File paired 'Bed))]
getSortedBed input = concatMap split $ concatMap split $
    input & mapped.replicates.mapped.files %~ map f . filter filterFn . lefts
  where
    filterFn x = getFileType x == Bed && x `hasTag` NameSorted && x `hasTag` Gzip
    f fl | fl `hasTag` PairedEnd = Right $ fromSomeFile fl
         | otherwise = Left $ fromSomeFile fl
      
demultiplex :: SCATACSeqConfig config
            => SCATACSeq S ( (File '[Gzip] 'Fastq, File '[Gzip] 'Fastq)
                           , File '[Gzip] 'Fastq )
            -> ReaderT config IO
                (SCATACSeq S (File '[Demultiplexed, Gzip] 'Fastq, File '[Demultiplexed, Gzip] 'Fastq))
demultiplex input = do
    bcLen <- fromMaybe (error "cell barcode length was not provided") <$> asks _scatacseq_cell_barcode_length
    dir <- asks ((<> "/Fastq") . _scatacseq_output_dir) >>= getPath
    let output1 = printf "%s/%s_rep%d_demulti_R1.fastq.gz" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
        output2 = printf "%s/%s_rep%d_demulti_R2.fastq.gz" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ ( \((fq1, fq2), fqidx) -> liftIO $ withTemp Nothing $ \tmp -> do
        stat <- runResourceT $ runConduit $ streamFastqGzip (fqidx^.location) .| takeC 100000000 .| barcodeStat bcLen
        B.writeFile tmp $ B.unlines $ map (B.pack . show . snd) $ U.toList stat
        thres <- fmap (read . T.unpack . head . T.lines) $ shelly $ run "taiji-utils" ["barcode", T.pack tmp]
        let bcMap = mkBarcodeMap $ map fst $ take (truncate thres) $ U.toList stat
        demulti output1 output2 bcMap bcLen (fqidx^.location) (fq1^.location) $ fq2^.location
        return ( location .~ output1 $ emptyFile
               , location .~ output2 $ emptyFile )
        )