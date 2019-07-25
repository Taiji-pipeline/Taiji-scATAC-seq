{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE QuasiQuotes       #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE DataKinds #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Diff where

import qualified Data.Vector as V
import qualified Data.ByteString.Char8 as B
import qualified Data.Text as T
import Bio.Data.Bed
import Data.Conduit.Internal (zipSources)
import System.Random.MWC (create)
import System.Random.MWC.Distributions (uniformShuffle)
import qualified Data.HashSet as S
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
          => (SCATACSeq S (File tags 'Other), FilePath)
          -> ReaderT config IO (SCATACSeq S (File '[] 'Tsv))
diffPeaks (input, ref) = do
    dir <- asks ((<> "/Diff/") . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_pvalues.txt" dir
            (T.unpack $ input^.eid) (input^.replicates._1)
    input & replicates.traversed.files %%~ liftIO . ( \fl -> do
        shelly $ run_ "sc_utils" [ "diff", T.pack output
            , "--fg", T.pack $ fl^.location 
            , "--bg", T.pack ref ]
        return $ location .~ output $ emptyFile )

getDiffPeaks :: SCATACSeqConfig config
             => ( File '[Gzip] 'NarrowPeak
                , SCATACSeq S (File '[] 'Tsv) )
             -> ReaderT config IO (SCATACSeq S (File '[Gzip] 'NarrowPeak))
getDiffPeaks (pk, input) = do
    dir <- asks ((<> "/Diff/Peak/") . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_diff.narrowpeak.gz" dir
            (T.unpack $ input^.eid) (input^.replicates._1)
    input & replicates.traversed.files %%~ liftIO . ( \fl -> do
        idx <- fmap S.fromList $ filterFDR $ fl^.location
        runResourceT $ runConduit $
            zipSources (iterateC succ 0) (streamBedGzip $ pk^.location :: ConduitT i NarrowPeak (ResourceT IO) ()) .|
            filterC ((`S.member` idx) . fst) .| mapC snd .| sinkFileBedGzip output
        return $ location .~ output $ emptyFile )

filterFDR :: FilePath -> IO [Int]
filterFDR fl = do
    (idx, prob) <- unzip . map (f . B.words) . B.lines <$> B.readFile fl
    prob' <- p_adjust prob
    return $ fst $ unzip $ filter ((<0.01) . snd) $ zip idx prob'
  where
    f [a,b] = (readInt a, readDouble b)

p_adjust :: [Double] -> IO [Double]
p_adjust xs = R.runRegion $ R.fromSomeSEXP <$>
    [r| p.adjust(xs_hs, method = "fdr") |]