{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE QuasiQuotes #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Diff
    ( specificPeaks
    ) where

import qualified Data.Vector as V
import qualified Data.ByteString.Char8 as B
import qualified Data.Text as T
import Bio.Data.Bed
import Bio.Utils.Functions (filterFDR)
import Shelly hiding (FilePath)

import Taiji.Prelude
import Taiji.Pipeline.SC.ATACSeq.Types
import qualified Taiji.Utils.DataFrame as DF

-- | Get cell-specific peaks
specificPeaks :: SCATACSeqConfig config
             => FilePath
             -> ( [(B.ByteString, File '[] 'NarrowPeak)]
                , (FilePath, FilePath, FilePath) )
             -> ReaderT config IO [(T.Text, File '[Gzip] 'NarrowPeak)]
specificPeaks prefix (peakFls, (_, scoreFl, pvalueFl)) = do
    dir <- asks ((<> asDir prefix) . _scatacseq_output_dir) >>= getPath
    liftIO $ do
        scores <- DF.readTable scoreFl
        pvalues <- DF.readTable pvalueFl
        let table = DF.zip scores pvalues
            regions = map mkPeak $ DF.rowNames table
        forM peakFls $ \(nm', peakFl) -> do
            let nm = T.pack $ B.unpack nm'
                values = V.toList $ table `DF.cindex` nm
                output = dir <> T.unpack nm <> "_diff_peaks.np.gz"
            peaks <- readBed $ peakFl^.location :: IO [BED3]
            candidates <- fmap (fst . unzip . V.toList . filterFDR fdr) $
                runConduit $ yieldMany (zipWith BEDExt regions values) .|
                    intersectBed peaks .| mapC f .| sinkVector
            runResourceT $ runConduit $ yieldMany candidates .| sinkFileBedGzip output 
            return (nm, location .~ output $ emptyFile)
  where
    f (BEDExt bed (sc, p)) = (npSignal .~ sc $ npPvalue .~ Just p $ convert bed, p)
    mkPeak :: T.Text -> BED3
    mkPeak x = let [chr, x'] = T.splitOn ":" x
                   [s,e] = T.splitOn "-" x'
               in asBed (B.pack $ T.unpack chr) (read $ T.unpack s) (read $ T.unpack e)
    fdr = 0.001