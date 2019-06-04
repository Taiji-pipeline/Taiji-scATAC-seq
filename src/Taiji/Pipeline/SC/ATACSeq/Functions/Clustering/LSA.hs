-- | Latent semantic analysis (LSA) is a technique in natural language
-- processing, in particular distributional semantics, of analyzing
-- relationships between a set of documents and the terms they contain
-- by producing a set of concepts related to the documents and terms.
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE DataKinds #-}

module Taiji.Pipeline.SC.ATACSeq.Functions.Clustering.LSA
    ( clust
    , performLSA
    ) where

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

{-
-- | Using LSA + graphical clustering
lsaClust :: (Elem 'Gzip tags ~ 'True, SCATACSeqConfig config)
         => SCATACSeq S (File tags 'Other)
         -> ReaderT config IO (SCATACSeq S [CellCluster])
lsaClust input = do
    dir <- asks ((<> "/ReadCount") . _scatacseq_output_dir) >>= getPath
    tmp <- asks _scatacseq_temp_dir
    let lsaNpy = printf "%s/%s_rep%d_readcount_lsa.npy" dir
            (T.unpack $ input^.eid) (input^.replicates._1)
    input & replicates.traversed.files %%~ liftIO . lsaClust' tmp lsaNpy

-- | Using LSA + graphical clustering
lsaClustMerged :: (Elem 'Gzip tags ~ 'True, SCATACSeqConfig config)
               => File tags 'Other
               -> ReaderT config IO [CellCluster]
lsaClustMerged input = do
    dir <- asks ((<> "/ReadCount") . _scatacseq_output_dir) >>= getPath
    tmp <- asks _scatacseq_temp_dir
    let lsaNpy = dir <> "/merged_lsa.npy"
    liftIO $ lsaClust' tmp lsaNpy input
    -}

    {-
-- | Using LSA + graphical clustering
lsaClust' :: Elem 'Gzip tags ~ 'True
          => Maybe FilePath    -- ^ temp dir
          -> FilePath    -- ^ LSA output
          -> File tags 'Other
          -> IO [CellCluster]
lsaClust' dir lsaNpy input = withTempDir dir $ \tmpD -> do
      tfidf (input^.location) $ tmpD <> "/tfidf"
      shelly $ run_ "sc_utils" ["run", T.pack tmpD <> "/tfidf", T.pack lsaNpy]
      shelly $ run_ "sc_utils" ["embed", T.pack lsaNpy, T.pack tmpD <> "/embed"]
      shelly $ run_ "sc_utils" ["clust", T.pack lsaNpy, T.pack tmpD <> "/clust"]

      clusters <- readClusters $ tmpD <> "/clust"
      sp <- mkSpMatrix readInt $ input^.location
      cellBcs <- runResourceT $ runConduit $
          streamRows sp .| mapC fst .| sinkList
      cells <- readCells cellBcs $ tmpD <> "/embed"
      return $ zipWith (\i -> CellCluster $ B.pack $ "C" ++ show i) [1::Int ..] $
          map (map (cells V.!)) clusters
  where
    readClusters fl = map (map readInt . B.split ',') . B.lines <$>
        B.readFile fl
    readCells bc fl = V.fromList . zipWith f bc .
        map (map readDouble . B.split '\t') . B.lines <$> B.readFile fl
      where
        f i [x,y,z] = Cell i x y z
        f _ _ = error "formatting error"
{-# INLINE lsaClust' #-}
-}

clust :: Maybe FilePath      -- ^ temp dir
      -> (File '[] 'Tsv, File '[Gzip] 'Tsv)   -- ^ lsa input matrix
      -> IO [CellCluster]
clust dir (coverage, mat) = withTempDir dir $ \tmpD -> do
      runResourceT $ runConduit $ seqDepthC .| mapC snd .|
          unlinesAsciiC .| sinkFile (tmpD <> "/coverage")
      shelly $ run_ "sc_utils" ["embed", T.pack $ mat^.location, T.pack tmpD <> "/embed"]
      shelly $ run_ "sc_utils" [ "clust", "--coverage"
          , T.pack tmpD <> "/coverage", T.pack $ mat^.location, T.pack tmpD <> "/clust" ]

      let sourceCells = getZipSource $ (,,) <$>
              ZipSource (iterateC succ 0) <*>
              ZipSource seqDepthC <*>
              ZipSource ( sourceFile (tmpD <> "/embed") .| linesUnboundedAsciiC .|
                mapC (map readDouble . B.split '\t') )
      cells <- runResourceT $ runConduit $ sourceCells .| mapC f .| sinkVector
      clusters <- readClusters $ tmpD <> "/clust"
      return $ zipWith (\i -> CellCluster $ B.pack $ "C" ++ show i) [1::Int ..] $
          map (map (cells V.!)) clusters
  where
    readClusters fl = map (map readInt . B.split ',') . B.lines <$>
        B.readFile fl
    seqDepthC = sourceFile (coverage^.location) .| linesUnboundedAsciiC .|
        mapC ((\[a,b] -> (a,b)) . B.split '\t')
    f (i, (bc, dep), [x,y,z]) = Cell i x y z bc $ readInt dep
    f _ = error "formatting error"
{-# INLINE clust #-}

performLSA :: (Elem 'Gzip tags ~ 'True, SCATACSeqConfig config)
           => SCATACSeq S (File tags 'Other)
           -> ReaderT config IO (SCATACSeq S (File '[] 'Tsv, File '[Gzip] 'Tsv))
performLSA input = do
    dir <- asks ((<> "/LSA") . _scatacseq_output_dir) >>= getPath
    tmp <- asks _scatacseq_temp_dir
    let output = printf "%s/%s_rep%d_lsa.tsv.gz" dir
            (T.unpack $ input^.eid) (input^.replicates._1)
        rownames = printf "%s/%s_rep%d_lsa.rownames.txt" dir
            (T.unpack $ input^.eid) (input^.replicates._1)
    input & replicates.traversed.files %%~ liftIO . ( \fl -> do
        lsa tmp output fl
        sp <- mkSpMatrix readInt $ fl^.location
        runResourceT $ runConduit $
            streamRows sp .| mapC f .| unlinesAsciiC .| sinkFile rownames
        return ( location .~ rownames $ emptyFile
               , location .~ output $ emptyFile ) )
  where
    f (nm, xs) = nm <> "\t" <> fromJust (packDecimal $ foldl1' (+) $ map snd xs)

-- | Perform LSA.
lsa :: Elem 'Gzip tags ~ 'True
    => Maybe FilePath    -- ^ temp dir
    -> FilePath    -- ^ LSA output
    -> File tags 'Other
    -> IO ()
lsa dir lsaNpy input = withTemp dir $ \tmp -> do
      tfidf (input^.location) tmp
      shelly $ run_ "sc_utils" ["svd", T.pack tmp, T.pack lsaNpy]
{-# INLINE lsa #-}

--------------------------------------------------------------------------------
-- Helper functions
--------------------------------------------------------------------------------

tfidf :: FilePath -> FilePath -> IO ()
tfidf input output = do
    sp <- mkSpMatrix readInt input
    idf <- mkIDF sp
    let header = B.pack $
            printf "Sparse matrix: %d x %d" (_num_row sp) (_num_col sp)
    runResourceT $ runConduit $
        streamRows sp .| mapC (second (\x -> idf `dotProd` binarizeTF x)) .|
        (yield header >> mapC showRow) .| unlinesAsciiC .| gzip .| sinkFile output
{-# INLINE tfidf #-}

-- | The weighting functions transform each cell, a_{ij} of A, to be the
-- product of a local term weight, l_{ij}, which describes the relative
-- frequency of a term in a document, and a global weight, g_{i}, which
-- describes the relative frequency of the term within the entire
-- collection of documents.
type IDF = U.Vector Double

mkIDF :: SpMatrix a -> IO IDF
mkIDF sp = do
    tc <- runResourceT $ runConduit $ streamRows sp .| countTerms
    return $ U.map
        (\x -> log (fromIntegral nDoc + 1) - log (fromIntegral x + 1)) tc
  where
    nDoc = _num_row sp
    countTerms = do
        vec <- lift $ UM.replicate (_num_col sp) (0 :: Int)
        mapM_C $ mapM_ (UM.unsafeModify vec (+1) . fst) . snd
        U.unsafeFreeze vec

binarizeTF :: [(Int, Int)] -> [(Int, Double)]
binarizeTF = map $ second (\x -> if x > 0 then 1 else 0)

normalizeTF :: [(Int, Double)] -> [(Int, Double)]
normalizeTF xs = map (second (\x -> log x - n)) xs
  where
    n = log $ (foldl1' (+) $ map snd xs) / 10000

showRow :: Row Double -> B.ByteString
showRow (name, xs) = B.intercalate "\t" $ name : map f xs
  where
    f (i,v) = fromJust (packDecimal i) <> "," <> toShortest v

-- | Dot product between vectors
dotProd :: U.Vector Double -> [(Int, Double)] -> [(Int, Double)]
dotProd vec = map $ \(i, x) -> (i, x * vec U.! i)
{-# INLINE dotProd #-}