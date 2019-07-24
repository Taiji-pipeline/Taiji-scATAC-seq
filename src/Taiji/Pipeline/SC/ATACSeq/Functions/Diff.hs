{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE DataKinds #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Diff where

import qualified Data.Vector as V
import qualified Data.ByteString.Char8 as B
import qualified Data.Text as T
import System.Random.MWC (create)
import System.Random.MWC.Distributions (uniformShuffle)
import qualified Data.HashSet as S
import Data.Conduit.Zlib (gzip)
import Shelly hiding (FilePath)

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