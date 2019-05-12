{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE OverloadedStrings #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Utils
    ( CellCluster(..)
    , readCellCluster
    , writeCellCluster
    , extractBarcode
    , readPromoters
    , getGenomeIndex
    , CutSite(..)
    , CutSiteIndex
    , withCutSiteIndex
    , getKeys
    , lookupIndex
    , createCutSiteIndex
    ) where

import Conduit
import Bio.Data.Bed
import System.IO
import Data.Conduit.Internal (zipSinks)
import           Bio.RealWorld.GENCODE
import Control.Arrow (second)
import Control.Lens
import Data.Either (either)
import qualified Data.Map.Strict as M
import qualified Data.HashMap.Lazy as HM
import qualified Data.ByteString.Char8 as B
import qualified Data.ByteString as BS
import Data.Serialize (encode, decode)
import System.IO.Temp (withTempFile)
import Control.DeepSeq (force)
import Control.Exception (bracket)
import           Data.List.Ordered       (nubSort)
import Control.Monad.Reader (asks)
import           Bio.Seq.IO
import Data.Maybe
import Control.Monad (unless)
import           System.FilePath               (takeDirectory)
import           Shelly                        (fromText, mkdir_p, shelly,
                                                test_f)
import Scientific.Workflow
import Data.CaseInsensitive (CI)
import qualified Data.Text as T

import Taiji.Pipeline.SC.ATACSeq.Types

data CellCluster = CellCluster
    { _cluster_name :: B.ByteString
    , _cluster_member :: [B.ByteString] }

writeCellCluster :: FilePath -> [CellCluster] -> IO ()
writeCellCluster output = B.writeFile output . B.unlines . map showCellCluster

readCellCluster :: FilePath -> IO [CellCluster]
readCellCluster input = map f . B.lines <$> B.readFile input
  where
    f xs = let [nm, members] = B.split '\t' xs
           in CellCluster nm $ B.split ',' members

showCellCluster :: CellCluster -> B.ByteString
showCellCluster CellCluster{..} = _cluster_name <> "\t" <>
    B.intercalate "," _cluster_member

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

getGenomeIndex :: SCATACSeqConfig config => WorkflowConfig config FilePath
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
    n <- either error id . decode <$> BS.hGet h 8
    (locationMap, chrs) <- either error id . decode <$> BS.hGet h n
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
encodeFile output chrs = mapC (second (encodeCutSite chrs . map bedToCutSite)) .|
    fmap snd (zipSinks sink1 sink2)
  where
    sink1 = mapC snd .| sinkFile output
    sink2 = snd <$> foldlC f (0, M.empty)
      where
        f (pos, acc) (nm, bs) = force (pos + fromIntegral n, M.insert nm (pos, n) acc)
          where 
            n = BS.length bs
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