{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE BangPatterns #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Utils
    ( extractBarcode
    , readPromoters
    , CutSite(..)
    , CutSiteIndex
    , withCutSiteIndex
    , getKeys
    , lookupIndex
    , createCutSiteIndex

      -- * Sparse Matrix
    , SpMatrix(..)
    , Row
    , mkSpMatrix
    , streamRows
    , sinkRows
    , decodeRowWith
    , encodeRowWith

    , filterCols
    , concatMatrix

      -- * Feature matrix
    , mkFeatMat
    , groupCells
    , mergeMatrix

    , computeSS
    , computeRAS
    , computeCDF
    ) where

import Bio.Data.Bed
import Bio.Data.Bed.Types
import System.IO
import Data.Conduit.List (groupBy)
import Data.Conduit.Internal (zipSinks)
import           Bio.RealWorld.GENCODE
import Control.Arrow (first, second)
import Data.Either (either)
import qualified Data.Map.Strict as M
import qualified Data.HashMap.Strict as HM
import qualified Data.Set as S
import qualified Data.ByteString.Char8 as B
import qualified Data.ByteString.Lazy as BS
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector as V
import qualified Data.Vector.Mutable as VM
import Data.ByteString.Lex.Integral (packDecimal)
import Bio.Utils.Misc (readInt)
import Data.Conduit.Zlib (multiple, ungzip, gzip)
import Data.Binary (encode, decode)
import Data.Singletons.Prelude (Elem)
import System.IO.Temp (withTempFile)
import System.Random.MWC (create)
import System.Random.MWC.Distributions (normal)
import Control.DeepSeq (force)
import Control.Exception (bracket)
import qualified Data.Matrix as Mat
import           Data.List.Ordered       (nubSort)
import Data.CaseInsensitive (CI)
import qualified Data.Text as T
import Statistics.Sample (varianceUnbiased, mean)

import Taiji.Prelude hiding (groupBy)
import qualified Taiji.Utils.DataFrame as DF

-------------------------------------------------------------------------------
-- BASIC
-------------------------------------------------------------------------------

-- | Extract barcode information from sequence name
extractBarcode :: B.ByteString -> B.ByteString
extractBarcode = head . B.split ':'
{-# INLINE extractBarcode #-}

-- | Get a list of potential TSS from GTF file
readPromoters :: FilePath -> IO (BEDTree [CI B.ByteString])
readPromoters = fmap (bedToTree (++) . concatMap fn) . readGenes
  where
    fn :: Gene -> [(BED3, [CI B.ByteString])]
    fn Gene{..} = map g $ nubSort tss
      where
        g x | geneStrand = (asBed geneChrom (max 0 $ x - 5000) (x + 1000), [geneName])
            | otherwise = (asBed geneChrom (max 0 $ x - 1000) (x + 5000), [geneName])
        tss | geneStrand = geneLeft : map transLeft geneTranscripts
            | otherwise = geneRight : map transRight geneTranscripts
{-# INLINE readPromoters #-}

-------------------------------------------------------------------------------
-- CutSiteIndex
-------------------------------------------------------------------------------

data CutSiteIndex = CutSiteIndex
    { _file_header :: CutSiteIndexHeader
    , _file_handle :: Handle }

data CutSiteIndexHeader = CutSiteIndexHeader
    { _header_location_map :: LocationMap
    , _chroms :: [B.ByteString]
    , _off_set :: Int }

data CutSite = CutSite
    { _cut_site_chr :: B.ByteString
    , _cut_site_pos :: Int }

type LocationMap = M.Map B.ByteString (Integer, Int)

getKeys :: CutSiteIndex -> [B.ByteString]
getKeys = M.keys . _header_location_map . _file_header
{-# INLINE getKeys #-}

lookupIndex :: B.ByteString -> CutSiteIndex -> IO (Maybe [CutSite])
lookupIndex key idx = case lookupHeader key (_file_header idx) of
    Nothing -> return Nothing
    Just (pos, n) -> do
        hSeek (_file_handle idx) AbsoluteSeek pos
        Just . decodeCutSite (_chroms $ _file_header idx) <$>
            BS.hGet (_file_handle idx) n
{-# INLINE lookupIndex #-}

lookupHeader :: B.ByteString
             -> CutSiteIndexHeader -> Maybe (Integer, Int)
lookupHeader key header = do
    (pos, n) <- M.lookup key $ _header_location_map header
    return (pos + fromIntegral (_off_set header), n)
{-# INLINE lookupHeader #-}

withCutSiteIndex :: FilePath -> (CutSiteIndex -> IO a) -> IO a
withCutSiteIndex idx = bracket (openCutSiteIndex idx) closeCutSiteIndex
{-# INLINE withCutSiteIndex #-}

openCutSiteIndex :: FilePath -> IO CutSiteIndex
openCutSiteIndex idx = do
    h <- openFile idx ReadMode
    n <- decode <$> BS.hGet h 8
    (locationMap, chrs) <- decode <$> BS.hGet h n
    return $ CutSiteIndex (CutSiteIndexHeader locationMap chrs (n + 8)) h
{-# INLINE openCutSiteIndex #-}

closeCutSiteIndex :: CutSiteIndex -> IO ()
closeCutSiteIndex = hClose . _file_handle
{-# INLINE closeCutSiteIndex #-}

createCutSiteIndex :: FilePath
                   -> [B.ByteString]
                   -> ConduitT () (B.ByteString, [BED]) (ResourceT IO) ()
                   -> IO ()
createCutSiteIndex output chrs source = withTempFile "./" "tmp_index" $ \tmp h_tmp -> do
    hClose h_tmp
    locationMap <- runResourceT $ runConduit $ source .| encodeFile tmp chrs
    let header = encode (locationMap, chrs)
    withFile output WriteMode $ \h -> do
        BS.hPutStr h $ encode (BS.length header) <> header
        runResourceT $ runConduit $ sourceFileBS tmp .| sinkHandle h
{-# INLINE createCutSiteIndex #-}

encodeFile :: (MonadResource m, PrimMonad m) 
           => FilePath
           -> [B.ByteString]
           -> ConduitT (B.ByteString, [BED]) Void m LocationMap
encodeFile output chrs = mapC (second (BS.toStrict . encodeCutSite chrs . map bedToCutSite)) .|
    fmap snd (zipSinks sink1 sink2)
  where
    sink1 = mapC snd .| sinkFile output
    sink2 = snd <$> foldlC f (0, M.empty)
      where
        f (pos, acc) (nm, bs) = force (pos + fromIntegral n, M.insert nm (pos, n) acc)
          where 
            n = B.length bs
{-# INLINE encodeFile #-}

bedToCutSite :: BED -> CutSite
bedToCutSite bed
    | bed^.strand == Just False = CutSite (bed^.chrom) (bed^.chromEnd)
    | otherwise = CutSite (bed^.chrom) (bed^.chromStart)
{-# INLINE bedToCutSite #-}

encodeCutSite :: [B.ByteString] -> [CutSite] -> BS.ByteString
encodeCutSite chrs sites = encode $ flip map chrs $ \chr ->
    map _cut_site_pos $ HM.lookupDefault [] chr $
    HM.fromListWith (++) $ map (\x -> (_cut_site_chr x, [x])) sites
{-# INLINE encodeCutSite #-}

decodeCutSite :: [B.ByteString] -> BS.ByteString -> [CutSite]
decodeCutSite chrs = concatMap f . zip chrs . either error id . decode
  where
    f (chr, xs) = map (\x -> CutSite chr x) xs
{-# INLINE decodeCutSite #-}

-------------------------------------------------------------------------------
-- Sparse Matrix
-------------------------------------------------------------------------------

data SpMatrix a = SpMatrix
    { _num_row :: Int
    , _num_col :: Int
    , _filepath :: FilePath
    , _decoder :: FilePath -> ConduitT () (Row a) (ResourceT IO) ()
    }

type Row a = (B.ByteString, [(Int, a)])

mkSpMatrix :: (B.ByteString -> a)   -- ^ Element decoder
           -> FilePath -> IO (SpMatrix a)
mkSpMatrix f input = do
    header <- runResourceT $ runConduit $ sourceFile input .| multiple ungzip .|
        linesUnboundedAsciiC .| headC
    case header of
        Nothing -> error "empty file"
        Just x -> do
            let [n, m] = map (read . T.unpack . T.strip) $ T.splitOn "x" $
                    last $ T.splitOn ":" $ T.pack $ B.unpack x
            return $ SpMatrix n m input decodeSpMatrix
  where
    decodeSpMatrix x = sourceFile x .| multiple ungzip .|
        linesUnboundedAsciiC .| (headC >> mapC (decodeRowWith f))
{-# NOINLINE mkSpMatrix #-}

streamRows :: SpMatrix a -> ConduitT () (Row a) (ResourceT IO) ()
streamRows sp = (_decoder sp) (_filepath sp)

sinkRows :: Int   -- ^ Number of rows
         -> Int   -- ^ Number of cols
         -> (a -> B.ByteString) 
         -> FilePath
         -> ConduitT (Row a) Void (ResourceT IO) ()
sinkRows n m encoder output = do
    (l, _) <- (yield header >> mapC (encodeRowWith encoder)) .| zipSinks lengthC sink
    when (l /= n + 1) $ error "incorrect number of rows"
  where
    header = B.pack $ printf "Sparse matrix: %d x %d" n m
    sink = unlinesAsciiC .| gzip .| sinkFile output


decodeRowWith :: (B.ByteString -> a) -> B.ByteString -> Row a
decodeRowWith decoder x = (nm, map f values)
  where
    (nm:values) = B.split '\t' x
    f v = let [i, a] = B.split ',' v
          in (readInt i, decoder a)
{-# INLINE decodeRowWith #-}

encodeRowWith :: (a -> B.ByteString) -> Row a -> B.ByteString
encodeRowWith encoder (nm, xs) = B.intercalate "\t" $ nm : map f xs
  where
    f (i,v) = fromJust (packDecimal i) <> "," <> encoder v
{-# INLINE encodeRowWith #-}

filterCols :: FilePath   -- ^ New matrix
           -> [Int]      -- ^ Columns to be removed
           -> FilePath   -- ^ Old matrix
           -> IO ()
filterCols output idx input = do
    mat <- mkSpMatrix id input
    let header = B.pack $ printf "Sparse matrix: %d x %d" (_num_row mat) (_num_col mat - length idx)
        newIdx = U.fromList $ zipWith (-) [0 .. _num_col mat - 1] $ go (-1,0) (sort idx)
        f = map (first (newIdx U.!)) . filter (not . (`S.member` idx') . fst)
        idx' = S.fromList idx
    runResourceT $ runConduit $ streamRows mat .| mapC (second f) .|
        (yield header >> mapC (encodeRowWith id)) .| unlinesAsciiC .|
        gzip .| sinkFile output
  where
    go (prev, c) (i:x) = replicate (i-prev) c ++ go (i, c+1) x
    go (_, c) [] = repeat c

-- | Combine rows of matrices. 
concatMatrix :: FilePath   -- ^ Output merged matrix
             -> [(Maybe B.ByteString, FilePath)] -- ^ A list of matrix
             -> IO ()
concatMatrix output inputs = do
    mats <- forM inputs $ \(nm, fl) -> do
        mat <- mkSpMatrix id fl
        return (nm, mat)
    runResourceT $ runConduit $ merge mats .| sinkFile output
  where
    merge mats
        | any (/=nBin) (map (_num_col . snd) mats) = error "Column unmatched!"
        | otherwise = source .| (yield header >> mapC (encodeRowWith id)) .|
            unlinesAsciiC .| gzip
      where
        source = forM_ mats $ \(nm, mat) ->
            let f x = case nm of
                    Nothing -> x
                    Just n -> n <> "+" <> x
            in streamRows mat .| mapC (first f)
        nCell = foldl' (+) 0 $ map (_num_row . snd) mats
        nBin = _num_col $ snd $ head mats
        header = B.pack $ printf "Sparse matrix: %d x %d" nCell nBin

-- | Generating the cell by feature count matrix as a gzipped stream.
-- The input stream is raw tags grouped by cells.
mkFeatMat :: (PrimMonad m, MonadThrow m)
          => Int   -- ^ the number of cells
          -> [[BED3]]    -- ^ a list of regions
          -> ConduitT (B.ByteString, [BED]) B.ByteString m ()
mkFeatMat nCell regions = source .| unlinesAsciiC .| gzip
  where
    source = yield header >> mapC
        (encodeRowWith (fromJust . packDecimal) . second (countEachCell bedTree))
    bedTree = bedToTree (++) $ concat $
        zipWith (\xs i -> zip xs $ repeat [i]) regions [0::Int ..]
    nBin = length regions
    header = B.pack $ printf "Sparse matrix: %d x %d" nCell nBin
    countEachCell :: BEDTree [Int] -> [BED] -> [(Int, Int)]
    countEachCell beds = HM.toList . foldl' f HM.empty
      where
        f m bed = foldl' (\x k -> HM.insertWith (+) k (1::Int) x) m $
            concat $ intersecting beds query
          where
            query = case bed^.strand of
                Just False -> BED3 (bed^.chrom) (bed^.chromEnd - 1) (bed^.chromEnd)
                _ -> BED3 (bed^.chrom) (bed^.chromStart) (bed^.chromStart + 1)

groupCells :: Monad m => ConduitT BED (B.ByteString, [BED]) m ()
groupCells = groupBy ((==) `on` (^.name)) .|
    mapC (\xs -> (fromJust $ head xs^.name, xs))
{-# INLINE groupCells #-}

-- | Merge sparse matrix.
mergeMatrix :: Elem 'Gzip tags1 ~ 'True
            => [(B.ByteString, (File tags1 file, File tags2 'Other))]  -- ^ (regions, matrix)
            -> FilePath   -- ^ Index file output
            -> ConduitT () B.ByteString (ResourceT IO) ()
mergeMatrix inputs idxOut = do
    indices <- liftIO $ getIndices $ map (fst . snd) inputs
    liftIO $ runResourceT $ runConduit $ yieldMany indices .| sinkFileBedGzip idxOut 
    let idxMap = M.fromList $ zip indices [0..]
    mats <- liftIO $ forM inputs $ \(nm, (idxFl, matFl)) -> do
        idxVec <- fmap (V.map (flip (M.findWithDefault undefined) idxMap)) $
            runResourceT $ runConduit $ streamBedGzip (idxFl^.location) .| sinkVector
        mat <- mkSpMatrix readInt $ matFl^.location
        return ( nm, mat
            { _decoder = \x -> _decoder mat x .|
                mapC (second (map (first (idxVec V.!))))
            , _num_col = M.size idxMap } )
    merge mats
  where
    merge ms
        | any (/=nBin) (map (_num_col . snd) ms) = error "Column unmatched!"
        | otherwise = source .|
            (yield header >> mapC (encodeRowWith (fromJust . packDecimal))) .|
            unlinesAsciiC .| gzip
      where
        source = forM_ ms $ \(sample, m) -> streamRows m .| 
            mapC (first (\x -> sample <> "+" <> x))
        nCell = foldl' (+) 0 $ map (_num_row . snd) ms
        nBin = _num_col $ snd $ head ms
        header = B.pack $ printf "Sparse matrix: %d x %d" nCell nBin
    getIndices :: [File tags file] -> IO [BED3]
    getIndices idxFls = fmap S.toList $ runResourceT $ runConduit $
        mapM_ (streamBedGzip . (^.location)) idxFls .|
        foldlC (flip S.insert) S.empty

computeRAS :: FilePath -> IO (V.Vector Double)
computeRAS fl = do
    mat <- mkSpMatrix readInt fl
    fmap accScore $ runResourceT $ runConduit $
        streamRows mat .| sink (_num_col mat)
  where
    sink n = do
        v <- lift $ VM.replicate n 0.1
        mapM_C $ \(_, xs) -> forM_ xs $ \(i, x) ->
            VM.unsafeModify v (+fromIntegral x) i
        lift $ V.unsafeFreeze v
    accScore xs = V.map (\x -> x * 1000000 / s) xs
      where
        s = V.sum xs

-- | Compute Specificity Score (SS).
computeSS :: DF.DataFrame Double -> DF.DataFrame Double
computeSS df = df{ DF._dataframe_data = mat}
  where
    mat = Mat.fromRows $ map (q_t . normalize) $ Mat.toRows $
        DF._dataframe_data df
    q_t xs = V.map (\x -> e - logBase 2 x) xs
      where
        e = negate $ V.sum $ V.map (\p -> p * logBase 2 p) xs
    normalize xs = V.map (/s) xs'
      where
        s = V.sum xs'
        xs' = V.map (+pseudoCount) xs
        pseudoCount = 1
{-# INLINE computeSS #-}

{-
buildNullModel :: DF.DataFrame Double -> IO (V.Vector Double)
buildNullModel df = do
    fmap (V.concat . map (q_t . normalize) . Mat.toRows) . shuffleMatrix
  where
    mat = let k = truncate $ 0.05 * fromIntegral (Mat.rows $ DF._dataframe_data df)
          in Mat.fromColumns $
              map (V.take k . V.modify (\x -> I.partialSortBy (flip compare) x k)) $
              Mat.toColumns $ DF._dataframe_data df
    shuffleMatrix :: Mat.Matrix Double -> IO (Mat.Matrix Double)
    shuffleMatrix mat = fmap (Mat.fromVector (Mat.dim mat)) $
        create >>= uniformShuffle (Mat.flatten mat)
    mkVector vec = 
      where
        vec' = V.modify (\x -> I.partialSortBy (flip compare) x k) vec
        k = truncate $ 0.1 * fromIntegral (V.length vec)
-}

computeCDF :: DF.DataFrame Double -> IO (V.Vector Double, Double, Double)
computeCDF df = do
    gen <- create
    let std = let rows = filter ((>5) . V.maximum) $ Mat.toRows $
                      Mat.map (+1) $ DF._dataframe_data df
                  f xs = (xs, 1 / entropy (normalize xs))
                  n = truncate $ (0.7 :: Double) * fromIntegral (length rows)
                  getFold xs = let m = mean xs in V.map (logBase 2 . (/m)) xs
                  entropy xs = negate $ V.sum $ V.map (\p -> p * logBase 2 p) xs
              in sqrt $ varianceUnbiased $ V.concat $ map (getFold . fst) $
                  take n $ sortBy (comparing snd) $ map f $ rows
    runConduit $ replicateMC nSample (V.toList <$> mkSample gen std) .|
        concatC .| mkCDF 0.001
  where
    ncol = Mat.cols $ DF._dataframe_data df
    nSample = 500000
    mkCDF res = do
        vec <- VM.replicate n 0
        mapM_C $ \x -> do
          let i = min n $ truncate $ x / res
          VM.unsafeModify vec (+1) i
        v <- V.scanl1' (+) <$> V.unsafeFreeze vec
        return (V.map (/(V.last v)) v, res, fromIntegral $ nSample*ncol)
      where
        n = truncate $ 2 * logBase 2 (fromIntegral ncol) / res
    -- | https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2396180/
    normalize xs = V.map (/s) xs
      where
        s = V.sum xs
    mkSample gen !std = do
        folds <- replicateM ncol $ normal 0 std gen 
        return $ q_t $ normalize $ V.fromList $ map (\x -> 2**x) folds
    q_t xs = V.map (\x -> e - logBase 2 x) xs
      where
        e = negate $ V.sum $ V.map (\p -> p * logBase 2 p) xs
{-# INLINE computeCDF #-}