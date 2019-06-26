{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE DataKinds #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Clustering.LDA (performLDA) where 

import Control.Arrow
import Data.Conduit.Zlib (gzip)
import Data.ByteString.Lex.Integral (packDecimal)
import Data.Double.Conversion.ByteString (toShortest)
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector as V
import Data.Singletons.Prelude (Elem)
import Bio.Utils.Misc (readDouble, readInt)
import qualified Data.ByteString.Char8 as B
import qualified Data.Vector.Unboxed.Mutable as UM
import qualified Data.Text as T
import Shelly (shelly, run_)

import qualified Language.R                        as R
import           Language.R.QQ
import qualified Data.Vector.SEXP as V
import Language.R.HExp

import Taiji.Prelude
import Taiji.Pipeline.SC.ATACSeq.Types
import Taiji.Pipeline.SC.ATACSeq.Functions.Utils

performLDA :: (Elem 'Gzip tags ~ 'True, SCATACSeqConfig config)
           => FilePath  -- ^ directory
           -> SCATACSeq S (File tags 'Other)
           -> ReaderT config IO (SCATACSeq S (File '[] 'Tsv, File '[Gzip] 'Tsv))
performLDA prefix input = do
    dir <- asks ((<> asDir prefix) . _scatacseq_output_dir) >>= getPath
    tmp <- asks _scatacseq_temp_dir
    let output = printf "%s/%s_rep%d_lda.tsv.gz" dir
            (T.unpack $ input^.eid) (input^.replicates._1)
        rownames = printf "%s/%s_rep%d_lda.rownames.txt" dir
            (T.unpack $ input^.eid) (input^.replicates._1)
    input & replicates.traversed.files %%~ liftIO . ( \fl -> do
        lda tmp output fl
        sp <- mkSpMatrix readInt $ fl^.location
        runResourceT $ runConduit $
            streamRows sp .| mapC f .| unlinesAsciiC .| sinkFile rownames
        return ( location .~ rownames $ emptyFile
               , location .~ output $ emptyFile ) )
  where
    f (nm, xs) = nm <> "\t" <> fromJust (packDecimal $ foldl1' (+) $ map snd xs)

-- | Perform LDA.
lda :: Elem 'Gzip tags ~ 'True
    => Maybe FilePath    -- ^ temp dir
    -> FilePath    -- ^ LDA output
    -> File tags 'Other
    -> IO ()
lda dir output input = withTemp dir $ \tmp -> do
      preprocess (input^.location) tmp
      shelly $ run_ "sc_utils" ["reduce", "--method", "lda", T.pack tmp, T.pack output]
{-# INLINE lda #-}

preprocess :: FilePath -> FilePath -> IO ()
preprocess input output = do
    sp <- mkSpMatrix readInt input
    let header = B.pack $ printf "Sparse matrix: %d x %d"
            (_num_row sp) (_num_col sp)
    runResourceT $ runConduit $
        streamRows sp .| mapC (second binarizeTF) .|
        (yield header >> mapC showRow) .| unlinesAsciiC .| gzip .| sinkFile output
{-# INLINE preprocess #-}

binarizeTF :: [(Int, Int)] -> [(Int, Double)]
binarizeTF = map $ second (\x -> if x > 0 then 1 else 0)

showRow :: Row Double -> B.ByteString
showRow (name, xs) = B.intercalate "\t" $ name : map f xs
  where
    f (i,v) = fromJust (packDecimal i) <> "," <> toShortest v


{-
-- | Run the LDA model implemented in cisTopic.
-- Input matrix: a precomputed matrix with cells as columns, regions as cells
-- and fragments/reads counts as values.
-- The rownames of these matrix must contain the region coordinates in
-- position format (e.g. chr1:123456-134567)
cisTopic :: FilePath    -- ^ input matrix
         -> IO ()
cisTopic input = R.runRegion $ do
    _ <- [r| library("cisTopic")
             cisTopicObject <- createcisTopicObject(counts_mel, project.name='cis')
             cisTopicObject <- runModels(cisTopicObject,
                topic=30, seed=987, nCores=13, burnin = 120, iterations = 150,
                addModels=FALSE)
            cisTopicObject <- selectModel(cisTopicObject)
-}
