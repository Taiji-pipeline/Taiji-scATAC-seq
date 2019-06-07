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

import Taiji.Prelude
import Taiji.Pipeline.SC.ATACSeq.Types
import Taiji.Pipeline.SC.ATACSeq.Functions.Utils

performLDA :: (Elem 'Gzip tags ~ 'True, SCATACSeqConfig config)
           => SCATACSeq S (File tags 'Other)
           -> ReaderT config IO (SCATACSeq S (File '[] 'Tsv, File '[Gzip] 'Tsv))
performLDA input = do
    dir <- asks ((<> "/LDA") . _scatacseq_output_dir) >>= getPath
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

