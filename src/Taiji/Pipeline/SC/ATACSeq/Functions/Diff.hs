{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE QuasiQuotes       #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE DataKinds #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Diff
    ( sampleCells
    , mkRefMat
    , diffPeaks
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
import Bio.Utils.Misc (readDouble, readInt)
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
         => ( [[B.ByteString]]
            , [SCATACSeq S (File '[Gzip] 'Other)] )
         -> ReaderT config IO FilePath
mkRefMat (bcs, [fl]) = do
    dir <- asks _scatacseq_output_dir >>= getPath . (<> "/Feature/Peak/")
    let output = dir <> "ref_cell_by_peak.mat.gz"
    liftIO $ do
        mat <- mkSpMatrix id $ fl^.replicates._2.files.location
        let header = B.pack $ printf "Sparse matrix: %d x %d" (S.size bc) (_num_col mat)
        runResourceT $ runConduit $ streamRows mat .|
            filterC ((`S.member` bc) . fst) .|
            (yield header >> mapC (encodeRowWith id)) .| unlinesAsciiC .|
            gzip .| sinkFile output
        return output
  where
    bc = S.fromList $ concat bcs
mkRefMat _ = undefined

diffPeaks :: SCATACSeqConfig config
          => ( SCATACSeq S (File tags 'Other)
             , File '[Gzip] 'NarrowPeak
             , FilePath )
          -> ReaderT config IO (SCATACSeq S (File '[] 'Tsv))
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

diffAnalysis :: FilePath  -- ^ foreground
             -> FilePath  -- ^ background
             -> IO [(Int, Double, Double, Double)]
diffAnalysis fg bk = withTemp Nothing $ \tmp -> do
    shelly $ run_ "sc_utils" [ "diff", T.pack tmp, "--fg", T.pack fg
        , "--bg", T.pack bk ]
    map (f . B.split '\t') . B.lines <$> B.readFile tmp
  where
    f [a,b,c,d] = (readInt a, readDouble b, readDouble c, readDouble d)
    f _ = error "wrong format!"