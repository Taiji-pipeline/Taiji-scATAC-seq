{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE DataKinds #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.ChromVar
    ( preChromVar
    , runChromVar
    ) where

import qualified Data.ByteString.Char8 as B
import Bio.Data.Bed
import Data.Binary (encodeFile, decodeFile)
import Bio.Data.Bed.Utils (fetchSeq)
import Data.Conduit.Internal (zipSinks)
import Data.Conduit.List (chunksOf)
import Bio.Seq hiding (length)
import Bio.Seq.IO (withGenome)
import qualified Data.Vector.Unboxed as U
import qualified Data.IntSet as IS
import System.Random.MWC
import Bio.ChromVAR

import Taiji.Prelude
import Taiji.Utils
import Taiji.Pipeline.SC.ATACSeq.Types

preChromVar :: SCATACSeqConfig config
            => ( Maybe FilePath
               , Maybe (File '[Gzip] 'NarrowPeak)
               , [SCATACSeq S (File '[Gzip] 'Other)] )
            -> ReaderT config IO
                [(Int, FilePath, File '[Gzip] 'Other, File '[] 'Other)]
preChromVar (Just motifFl, Just peakFl, inputs) = do
    dir <- asks _scatacseq_output_dir >>= getPath . (<> "/ChromVar/")
    genome <- asks ( fromMaybe (error "Genome index file was not specified!") .
        _scatacseq_genome_index )
    let f (Left _) = 0.5
        f (Right s) = gcContent (s :: DNA IUPAC)
        peakMat = dir <> "merged_cell_by_peak.mat.gz"
        bin = dir <> "chromVar.bin"
    liftIO $ do
        concatMatrix peakMat $ zip (repeat Nothing) $
            inputs^..folded.replicates._2.files.location
        mat <- mkSpMatrix readDouble peakMat
        gc <- withGenome genome $ \g -> runResourceT $ runConduit $ 
            streamBedGzip (peakFl^.location) .|
            mapC (convert :: NarrowPeak -> BED) .|
            fetchSeq g .| mapC f .| sinkVector
        readCount <- colSum mat
        let acc = U.map (logBase 10) readCount
        bk <- create >>= getBackgroundPeaks 50 (U.zip acc gc)
        encodeFile bin (readCount, bk)

        motifs <- decodeFile motifFl :: IO [(B.ByteString, IS.IntSet)]
        forM (zip (runIdentity $ runConduit $ yieldMany motifs .| chunksOf 10 .| sinkList) [0..]) $ \(ms, i) -> do
            let output = dir <> "motifs_" <> show i <> ".bin"
            encodeFile output ms
            return (i, output, location .~ peakMat $ emptyFile
                , location .~ bin $ emptyFile)
preChromVar _ = return []

runChromVar :: SCATACSeqConfig config
            => (Int, FilePath, File '[Gzip] 'Other, File '[] 'Other)
            -> ReaderT config IO ()
runChromVar (idx, motifFl, peakMat, bin) = do
    dir <- asks _scatacseq_output_dir >>= getPath . (<> "/ChromVar/")
    let output1 = dir <> show idx <> "_chromVar_deviation.tsv"
        output2 = dir <> show idx <> "_chromVar_zscore.tsv"
    liftIO $ do
        (names, motifs) <- unzip <$> decodeFile motifFl
        (readCount, pg) <- decodeFile $ bin^.location
        let expect = U.map (/(U.sum readCount)) readCount
        mat <- mkSpMatrix readDouble $ peakMat^.location
        let f xs = B.intercalate "\t" $ map toShortest xs
            header = B.intercalate "\t" names
            sink g output = (yield header >> (concatMapC g .| mapC f)) .|
                unlinesAsciiC .| sinkFile output
        _ <- runResourceT $ runConduit $ chromVar 1000 expect pg mat motifs .|
            zipSinks (sink fst output1) (sink snd output2)
        return ()
 
-- | The ChromVar algorithm.
chromVar :: Int   -- ^ Batch size
         -> U.Vector Double   -- ^ expectation
         -> [U.Vector Int]    -- ^ Background peak
         -> SpMatrix Double   -- ^ Cell by peak count matrix
         -> [IS.IntSet]   -- ^ motifs
         -> ConduitT () ([[Double]], [[Double]]) (ResourceT IO) () 
chromVar chunks expect bk peakMat motifs = do
    streamRows peakMat .| chunksOf chunks .| mapC (mkPeaks . snd . unzip) .|
        computeDeviation nMotif nPeak (U.toList expect) (map U.toList bk) motifs'
  where
    mkPeaks xs = (length xs, concat $ zipWith f [0..] xs)
      where
        f i x = map (\(j, v) -> (i,j,v)) x
    nMotif = length motifs
    nPeak = _num_col peakMat
    motifs' = concat $ zipWith f [0..] motifs
      where
        f j m = map (\i -> (i,j,1)) $ IS.toList m
{-# INLINE chromVar #-}