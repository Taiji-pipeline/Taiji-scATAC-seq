{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE OverloadedStrings #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Utils
    ( extractBarcode
    , readPromoters
    , getGenomeIndex
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
    , decodeRowWith
    , encodeRowWith

    , visualizeCluster
    ) where

import Bio.Data.Bed
import System.IO
import Data.Conduit.Internal (zipSinks)
import           Bio.RealWorld.GENCODE
import Control.Arrow (second)
import Data.Either (either)
import qualified Data.Map.Strict as M
import qualified Data.HashMap.Strict as HM
import qualified Data.ByteString.Char8 as B
import qualified Data.ByteString.Lazy as BS
import Data.ByteString.Lex.Integral (packDecimal)
import Bio.Utils.Misc (readInt)
import Data.Conduit.Zlib (multiple, ungzip)
import Data.Binary (encode, decode)
import System.IO.Temp (withTempFile)
import Control.DeepSeq (force)
import Control.Exception (bracket)
import           Data.List.Ordered       (nubSort)
import           Bio.Seq.IO
import           System.FilePath               (takeDirectory)
import           Shelly                        (fromText, mkdir_p, shelly,
                                                test_f)
import Data.CaseInsensitive (CI)
import qualified Data.Text as T

import Taiji.Prelude
import Taiji.Pipeline.SC.ATACSeq.Types
import Taiji.Utils.Plot
import Taiji.Utils.Plot.ECharts

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
        tss | geneStrand = geneLeft : map fst geneTranscripts
            | otherwise = geneRight : map snd geneTranscripts
{-# INLINE readPromoters #-}

getGenomeIndex :: SCATACSeqConfig config => ReaderT config IO FilePath
getGenomeIndex = do
    seqIndex <- asks ( fromMaybe (error "Genome index file was not specified!") .
        _scatacseq_genome_index )
    genome <- asks ( fromMaybe (error "Genome fasta file was not specified!") .
        _scatacseq_genome_fasta )
    shelly $ do
        fileExist <- test_f $ fromText $ T.pack seqIndex
        unless fileExist $ do
            mkdir_p $ fromText $ T.pack $ takeDirectory seqIndex
            liftIO $ mkIndex [genome] seqIndex
    return seqIndex
{-# INLINE getGenomeIndex #-}

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
            return $ SpMatrix n m input $ decodeSpMatrix f

streamRows :: SpMatrix a -> ConduitT () (Row a) (ResourceT IO) ()
streamRows sp = (_decoder sp) (_filepath sp)

decodeSpMatrix :: (B.ByteString -> a)
               -> FilePath -> ConduitT () (Row a) (ResourceT IO) ()
decodeSpMatrix f input = sourceFile input .| multiple ungzip .|
    linesUnboundedAsciiC .| (headC >> mapC (decodeRowWith f))
{-# INLINE decodeSpMatrix #-}

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

visualizeCluster :: FilePath
                 -> HM.HashMap B.ByteString Double -> [CellCluster] -> IO ()
visualizeCluster output marker cs = do
    savePlots output [] [scatter3D dat]
  where
    dat = flip map cs $ \(CellCluster nm cells) -> (B.unpack nm, map f cells) 
    f (Cell i x y z) = [x, y, z, HM.lookupDefault undefined i marker]