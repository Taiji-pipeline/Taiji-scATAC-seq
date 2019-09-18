{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE QuasiQuotes #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Diff
    ( sampleCells
    , mkRefMat
    , diffPeaks
    , specificPeaks
    ) where

import qualified Data.Vector as V
import qualified Data.ByteString.Char8 as B
import qualified Data.Text as T
import Bio.Data.Bed
import Bio.Utils.Functions (filterFDR)
import Data.Binary (decodeFile)
import Data.Conduit.Internal (zipSources)
import System.Random.MWC (create)
import System.Random.MWC.Distributions (uniformShuffle)
import qualified Data.HashSet as S
import qualified Data.HashMap.Strict as M
import Shelly hiding (FilePath)

import Taiji.Prelude
import Taiji.Pipeline.SC.ATACSeq.Types
import Taiji.Pipeline.SC.ATACSeq.Functions.Utils
import qualified Taiji.Utils.DataFrame as DF
import Taiji.Pipeline.SC.ATACSeq.Functions.Feature (streamMatrices)

sampleCells :: Int   -- ^ Number of cells
            -> SCATACSeq S (File tags 'Other)   -- ^ clusters
            -> IO (SCATACSeq S [[B.ByteString]])
sampleCells n input = input & replicates.traversed.files %%~ ( \fl -> do
    cls <- decodeFile $ fl^.location
    g <- create
    forM cls $ \cl -> do
        let v = V.fromList $ map _cell_barcode $ _cluster_member cl
        V.toList . V.take n <$> uniformShuffle v g
    )

mkRefMat :: SCATACSeqConfig config
         => FilePath   -- ^ Prefix
         -> [B.ByteString]   -- ^ barcodes
         -> [SCATACSeq S (File '[Gzip] 'Other)]    -- ^ matrices
         -> ReaderT config IO (Maybe (File '[Gzip] 'Other))
mkRefMat _ _ [] = return Nothing
mkRefMat prefix bcs mats = do
    dir <- asks ((<> asDir prefix) . _scatacseq_output_dir) >>= getPath
    let output = dir <> "/ref_mat.gz"
    liftIO $ do
        mat <- mkSpMatrix id $ head mats ^. replicates._2.files.location
        runResourceT $ runConduit $ streamMatrices id mats .|
            filterC ((`S.member` bcSet) . fst) .| sinkRows (S.size bcSet) (_num_col mat) id output
        return $ Just $ location .~ output $ emptyFile
  where
    bcSet = S.fromList bcs

diffPeaks :: SCATACSeqConfig config
          => FilePath 
          -> ( File '[Gzip] 'NarrowPeak
             , SCATACSeq S (File tags 'Other, File '[] 'NarrowPeak)   -- ^ Cluster peaks
             , File '[Gzip] 'Other ) -- ^ Ref matrix
          -> ReaderT config IO (SCATACSeq S (File '[Gzip] 'NarrowPeak))
diffPeaks prefix (peakFl, input, ref) = do
    dir <- asks ((<> asDir prefix) . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d.np.gz" dir
            (T.unpack $ input^.eid) (input^.replicates._1)
    input & replicates.traversed.files %%~ liftIO . ( \(fl, clPeak) -> do
        masterPeakList <- runResourceT $ runConduit $
            streamBedGzip (peakFl^.location) .| sinkList :: IO [BED3]
        peaks <- bedToTree const . map (\x -> (x, ())) <$>
            (readBed $ clPeak^.location :: IO [BED3])
        let idx = map fst $ filter (\(_, x) -> peaks `isIntersected` x) $
                zip [0..] masterPeakList
        stats <- fmap (M.fromList . map (\(a,b,c,d) -> (a, (b,c,d))) . filter (\x -> x^._4 < 0.01)) $
            diffAnalysis (fl^.location) (ref^.location) $ Just idx
        let f :: (Int, BED3) -> Maybe NarrowPeak
            f (i, peak) = case M.lookup i stats of
                Nothing -> Nothing
                Just (fold, pval, fdr) -> 
                    let logP | pval == 0 = 200
                             | otherwise = negate $ logBase 10 pval 
                        logFDR | fdr == 0 = 200
                               | otherwise = negate $ logBase 10 fdr
                    in Just $ npSignal .~ fold $
                        npPvalue .~ Just logP $
                        npQvalue .~ Just logFDR $ convert peak
        runResourceT $ runConduit $
            zipSources (iterateC succ 0) (streamBedGzip $ peakFl^.location) .|
            concatMapC f .| sinkFileBedGzip output
        return $ location .~ output $ emptyFile )


diffAnalysis :: FilePath  -- ^ foreground
             -> FilePath  -- ^ background
             -> Maybe [Int] -- ^ selected idx
             -> IO [(Int, Double, Double, Double)]
diffAnalysis fg bk idx = withTempDir Nothing $ \dir -> do
    let resFl = dir ++ "/res.txt"
        idxFl = dir ++ "/idx.txt"
        makeIdx x = do
            B.writeFile idxFl $ B.unlines $ map (B.pack . show) x
            return ["--index", T.pack idxFl]
    useIdx <- maybe (return []) makeIdx idx
    shelly $ run_ "sc_utils" $
        [ "diff", T.pack resFl, "--fg", T.pack fg
        , "--bg", T.pack bk ] ++ useIdx
    map (f . B.words) . B.lines <$> B.readFile resFl
  where
    f [a,b,c,d] = (readInt a, readDouble b, readDouble c, readDouble d)
    f _ = error "wrong format!"

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