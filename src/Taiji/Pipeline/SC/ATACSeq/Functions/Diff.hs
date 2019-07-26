{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE QuasiQuotes       #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE DataKinds #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Diff
    ( sampleCells
    , mkRefMat
    , diffPeaks
    , diffGenes
    ) where

import qualified Data.Vector as V
import qualified Data.ByteString.Char8 as B
import qualified Data.Text as T
import Bio.Data.Bed
import Data.Conduit.Internal (zipSources)
import System.Random.MWC (create)
import System.Random.MWC.Distributions (uniformShuffle)
import qualified Data.HashSet as S
import qualified Data.HashMap.Strict as M
import Data.Conduit.Zlib (gzip)
import Shelly hiding (FilePath)

import qualified Language.R                        as R
import           Language.R.QQ
import Language.R.HExp

import Taiji.Prelude
import Taiji.Pipeline.SC.ATACSeq.Types
import Taiji.Pipeline.SC.ATACSeq.Functions.Utils

sampleCells :: Int   -- ^ number of cells
            -> SCATACSeq S (File tags 'Other)
            -> IO [B.ByteString]
sampleCells n input = do
    mat <- mkSpMatrix id $ input^.replicates._2.files.location
    v <- runResourceT $ runConduit $ streamRows mat .| mapC fst .| sinkVector 
    g <- create
    map f . V.toList . V.take n <$> uniformShuffle v g
  where
    f x = B.pack (T.unpack $ input^.eid) <> "+" <> x

mkRefMat :: SCATACSeqConfig config
         => FilePath
         -> Bool
         -> ( [[B.ByteString]]
            , [SCATACSeq S (File '[Gzip] 'Other)] )
         -> ReaderT config IO FilePath
mkRefMat filename addName (bcs, fls) = do
    dir <- asks _scatacseq_output_dir >>= getPath
    let output = dir <> filename
    liftIO $ do
        mat <- mkSpMatrix id $ head fls^.replicates._2.files.location
        let header = B.pack $ printf "Sparse matrix: %d x %d" (S.size bc) (_num_col mat)
        runResourceT $ runConduit $ mapM_ sourceMat fls .|
            (yield header >> mapC (encodeRowWith id)) .| unlinesAsciiC .|
            gzip .| sinkFile output
        return output
  where
    bc = S.fromList $ concat bcs
    sourceMat input = do
        let filterFn x | addName = (B.pack (T.unpack $ input^.eid) <> "+" <> x) `S.member` bc
                       | otherwise = x `S.member` bc
        mat <- liftIO $ mkSpMatrix id $ input^.replicates._2.files.location
        streamRows mat .| filterC (filterFn . fst)

diffPeaks :: SCATACSeqConfig config
          => ( SCATACSeq S (File tags 'Other)
             , File '[Gzip] 'NarrowPeak
             , FilePath )
          -> ReaderT config IO (SCATACSeq S (File '[Gzip] 'NarrowPeak))
diffPeaks (input, peakFl, ref) = do
    dir <- asks ((<> "/Diff/Peak/") . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d.np.gz" dir
            (T.unpack $ input^.eid) (input^.replicates._1)
    input & replicates.traversed.files %%~ liftIO . ( \fl -> do
        stats <- fmap (M.fromList . map (\(a,b,c,d) -> (a, (b,c,d))) . filter (\x -> x^._4 < 0.01)) $
            diffAnalysis (fl^.location) ref
        let f (i, peak) = case M.lookup i stats of
                Nothing -> Nothing
                Just (fold, pval, fdr) -> Just $ npSignal .~ fold $
                    npPvalue .~ Just (negate $ logBase 10 pval) $
                    npQvalue .~ Just (negate $ logBase 10 fdr) $ peak
        runResourceT $ runConduit $
            zipSources (iterateC succ 0) (streamBedGzip $ peakFl^.location) .|
            concatMapC f .| sinkFileBedGzip output
        return $ location .~ output $ emptyFile )

{-
mkDiffPeakFig :: SCATACSeqConfig config 
              => [SCATACSeq S (File '[Gzip] 'NarrowPeak)]
              -> ReaderT config IO ()
mkDiffPeakFig inputs = do
    dir <- figDir
    forM_ inputs $ \input -> do
        input^.eid
        streamBedGzip (input^.replicates._2.files.location)
-}
    

diffGenes :: SCATACSeqConfig config
          => ( SCATACSeq S (File tags 'Other)
             , FilePath   -- ^ Gene names
             , FilePath ) -- ^ Ref matrix
          -> ReaderT config IO (SCATACSeq S (File '[] 'Tsv))
diffGenes (input, nameFl, ref) = do
    dir <- asks ((<> "/Diff/Gene/") . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d.tsv" dir
            (T.unpack $ input^.eid) (input^.replicates._1)
    input & replicates.traversed.files %%~ liftIO . ( \fl -> do
        stats <- fmap (M.fromList . map (\(a,b,c,d) -> (a, (b,c,d))) . filter (\x -> x^._4 < 0.01)) $
            diffAnalysis (fl^.location) ref
        let f (i, nm) = case M.lookup i stats of
                Nothing -> Nothing
                Just (fold, pval, fdr) -> Just $ B.intercalate "\t"
                    [ nm, toShortest fold, toShortest $ negate $ logBase 10 pval
                    , toShortest $ negate $ logBase 10 fdr ]
        runResourceT $ runConduit $ zipSources (iterateC succ 0)
            (sourceFile nameFl .| linesUnboundedAsciiC .| mapC (head . B.words)) .|
            concatMapC f .| unlinesAsciiC .| sinkFile output
        return $ location .~ output $ emptyFile )


diffAnalysis :: FilePath  -- ^ foreground
             -> FilePath  -- ^ background
             -> IO [(Int, Double, Double, Double)]
diffAnalysis fg bk = withTemp Nothing $ \tmp -> do
    shelly $ run_ "sc_utils" [ "diff", T.pack tmp, "--fg", T.pack fg
        , "--bg", T.pack bk ]
    map (f . B.words) . B.lines <$> B.readFile tmp
  where
    f [a,b,c,d] = (readInt a, readDouble b, readDouble c, readDouble d)
    f _ = error "wrong format!"

