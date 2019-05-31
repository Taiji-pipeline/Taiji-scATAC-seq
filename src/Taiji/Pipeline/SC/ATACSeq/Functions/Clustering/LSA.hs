-- | Latent semantic analysis (LSA) is a technique in natural language
-- processing, in particular distributional semantics, of analyzing
-- relationships between a set of documents and the terms they contain
-- by producing a set of concepts related to the documents and terms.
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE DataKinds #-}

module Taiji.Pipeline.SC.ATACSeq.Functions.Clustering.LSA
    (applyTFIDF) where

import Control.Monad.Primitive
import           Bio.Pipeline.Utils
import Control.Arrow
import Data.Conduit.Zlib (multiple, ungzip, gzip)
import Data.ByteString.Lex.Integral (packDecimal)
import Data.Double.Conversion.ByteString (toShortest)
import qualified Data.Vector.Unboxed as U
import Data.Singletons.Prelude (Elem)
import Bio.Utils.Misc (readDouble, readInt)
import qualified Data.ByteString.Char8 as B
import qualified Data.Vector.Unboxed.Mutable as UM
import qualified Data.Text as T
import Text.Printf (printf)

import Taiji.Prelude
import Taiji.Pipeline.SC.ATACSeq.Types

applyTFIDF :: (Elem 'Gzip tags ~ 'True, SCATACSeqConfig config)
           => SCATACSeq S (File tags 'Other)
           -> ReaderT config IO (SCATACSeq S (File tags 'Other))
applyTFIDF input = do
    dir <- asks ((<> "/ReadCount") . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_readcount_tfidf.txt.gz" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traversed.files %%~ liftIO . (\fl -> do
        tfidf (fl^.location) output
        return $ emptyFile & location .~ output )

tfidf :: FilePath -> FilePath -> IO ()
tfidf input output = do
    idf <- mkIDFModelFromFile input
    runResourceT $ runConduit $
        sourceFile input .| multiple ungzip .| linesUnboundedAsciiC .|
        transform idf .| unlinesAsciiC .| gzip .| sinkFile output
  where
    transform idf = headC >>= \case
        Nothing -> error "empty file"
        Just x -> do
            yield x
            mapC (showRow . second (applyIDF idf . binarizeTF) . getRow) 

-- | The weighting functions transform each cell, a_{ij} of A, to be the
-- product of a local term weight, l_{ij}, which describes the relative
-- frequency of a term in a document, and a global weight, g_{i}, which
-- describes the relative frequency of the term within the entire
-- collection of documents.
newtype IDFModel = IDFModel (U.Vector Double)

type Document a = [(Int, a)] 

binarizeTF :: Document Double -> Document Double
binarizeTF = map (second (\x -> if x > 0 then 1 else 0))

mkIDFModelFromFile :: FilePath -> IO IDFModel
mkIDFModelFromFile input = runResourceT $ runConduit $
    sourceFile input .| multiple ungzip .| linesUnboundedAsciiC .| sink
  where
    sink = headC >>= \case
        Nothing -> error "empty file"
        Just x -> do
            let n = readInt $ B.tail $ last $ B.split 'x' x
            mapC (snd . getRow) .| mkIDFModel n

getRow :: B.ByteString -> (B.ByteString, Document Double)
getRow x = (name, map f values)
  where
    (name:values) = B.split '\t' x
    f v = let [i, a] = B.split ',' v
          in (readInt i, readDouble a)

showRow :: (B.ByteString, Document Double) -> B.ByteString
showRow (name, xs) = B.intercalate "\t" $ name : map f xs
  where
    f (i,v) = fromJust (packDecimal i) <> "," <> toShortest v

mkIDFModel :: PrimMonad m => Int -> ConduitT (Document a) Void m IDFModel
mkIDFModel nTerm = do
    (nDoc, tc) <- getZipSink $ (,) <$> ZipSink lengthC <*> ZipSink countTerms
    return $ IDFModel $ U.map
        (\x -> log (fromIntegral (nDoc::Int) + 1) - log (fromIntegral x + 1)) tc
  where
    countTerms = do
        vec <- lift $ UM.replicate nTerm (0 :: Int)
        mapM_C $ \xs -> forM_ xs $ \(i,_) -> UM.unsafeModify vec (+1) i
        U.unsafeFreeze vec

applyIDF :: IDFModel -> Document Double -> Document Double
applyIDF (IDFModel vec) = map $ \(i, x) -> (i, x * vec U.! i)
{-# INLINE applyIDF #-}

{-
-- | The inverse document frequency is a measure of how much information
-- the word provides, i.e., if it's common or rare across all documents.
-- It is the logarithmically scaled inverse fraction of the documents that
-- contain the word (obtained by dividing the total number of documents by
-- the number of documents containing the term, and then taking the
-- logarithm of that quotient).
idf :: Int   -- ^ ith term
    -> TFMatrix -> Double
idf i (TFMatrix m) = logBase 2 $ nDoc / nDocWithTerm
  where
    nDocWithTerm = U.length $ U.filter (/=0) $ m `MU.takeRow` i
    nDoc = MU.cols m

tfidf :: Int -> Int -> TFMatrix -> Double
tfidf 
-}