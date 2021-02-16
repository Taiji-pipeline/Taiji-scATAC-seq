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

      -- * Feature matrix
    , mkFeatMat
    , groupCells
    ) where

import Bio.Data.Bed
import Bio.Data.Bed.Types
import System.IO
import Data.Conduit.List (groupBy)
import Data.Conduit.Internal (zipSinks)
import           Bio.RealWorld.GENCODE
import Control.Arrow (first, second)
import qualified Data.Map.Strict as M
import qualified Data.HashMap.Strict as HM
import qualified Data.Set as S
import qualified Data.ByteString.Char8 as B
import qualified Data.ByteString.Lazy as BS
import qualified Data.Vector as V
import Data.Conduit.Zlib (gzip)
import Data.Binary (encode, decode)
import Data.Singletons.Prelude (Elem)
import System.IO.Temp (withTempFile)
import Control.DeepSeq (force)
import Control.Exception (bracket)
import           Data.List.Ordered       (nubSort)
import Data.CaseInsensitive (CI)

import Taiji.Prelude hiding (groupBy)
import Taiji.Utils

-------------------------------------------------------------------------------
-- BASIC
-------------------------------------------------------------------------------

-- | Extract barcode information from sequence name
extractBarcode :: B.ByteString -> B.ByteString
extractBarcode = head . B.split ':'
{-# INLINE extractBarcode #-}

-- | Get a list of potential TSS from GTF file
readPromoters :: FilePath -> IO (BEDTree [CI B.ByteString])
readPromoters = fmap (bedToTree (++) . concatMap fn) . readGenesValidated
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