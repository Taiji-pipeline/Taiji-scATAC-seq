{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE GADTs #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Quantification
    ( mkCutSiteIndex
    , getBins
    , mkCellByBinMat
    , mergeMatrix
    , estimateExpr
    , mkExprTable
    ) where

import qualified Data.ByteString.Char8 as B
import Data.Conduit.List (groupBy)
import qualified Data.HashMap.Strict                  as M
import qualified Data.HashSet as S
import Data.ByteString.Lex.Integral (packDecimal)
import Data.Double.Conversion.ByteString (toShortest)
import           Bio.Utils.Misc (readInt, readDouble)
import Data.Conduit.Internal (zipSinks)
import Control.Arrow (first)
import Control.Monad.ST
import Bio.Data.Bed.Types
import Data.Conduit.Zlib (gzip)
import Bio.Data.Bed
import Bio.RealWorld.GENCODE (readGenes, Gene(..))
import           Data.CaseInsensitive  (mk, original, CI)
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Unboxed.Mutable as UM
import qualified Data.IntervalMap.Strict      as IM
import Data.Singletons.Prelude (Elem)
import           Bio.Pipeline.Utils
import Control.Arrow (second)
import           Data.List.Ordered                    (nubSort)
import qualified Data.Text as T
import Control.DeepSeq (force)
import Bio.Data.Bed.Utils
import Bio.Seq.IO (withGenome, getChrSizes)

import Taiji.Prelude hiding (groupBy)
import Taiji.Pipeline.SC.ATACSeq.Types
import Taiji.Pipeline.SC.ATACSeq.Functions.Utils

--------------------------------------------------------------------------------
-- Cut site map
--------------------------------------------------------------------------------

-- | Create the cut site map for every cell.
mkCutSiteIndex :: SCATACSeqConfig config
               => SCATACSeq S (File '[NameSorted, Gzip] 'Bed)
               -> ReaderT config IO (SCATACSeq S (File '[] 'Other))
mkCutSiteIndex input = do
    dir <- asks ((<> "/CutSiteIndex") . _scatacseq_output_dir) >>= getPath
    genome <- fromJust <$> asks _scatacseq_genome_index 
    let output = printf "%s/%s_rep%d.csidx" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . (\fl -> do
        chrs <- withGenome genome $ return . map fst . getChrSizes
        mkCutSiteIndex_ output fl chrs
        return $ emptyFile & location .~ output )
{-# INLINE mkCutSiteIndex #-}

mkCutSiteIndex_ :: FilePath
                -> File '[NameSorted, Gzip] 'Bed
                -> [B.ByteString]
                -> IO ()
mkCutSiteIndex_ output input chrs = createCutSiteIndex output chrs $
    streamBedGzip (input^.location) .| groupBy ((==) `on` (^.name)) .|
    mapC (\x -> (fromJust $ head x ^. name, x))
{-# INLINE mkCutSiteIndex_ #-}

--------------------------------------------------------------------------------
-- Cell by Bin Matrix
--------------------------------------------------------------------------------

-- | Get candidate bins.
getBins :: (Elem 'Gzip tags ~ 'True, SCATACSeqConfig config)
        => SCATACSeq S (File tags 'Bed)
        -> ReaderT config IO
            (SCATACSeq S (File tags 'Bed, File tags 'Bed, Int))
getBins input = do
    genome <- asks (fromJust . _scatacseq_genome_index)
    chrSize <- liftIO $ withGenome genome $ return . getChrSizes
    dir <- asks ((<> "/ReadCount") . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_bin_idx.bed.gz" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . (\fl -> do
        (n, rc) <- binCount fl chrSize res
        let lo = fromIntegral n * 0.001
            hi = fromIntegral n * 1
        runResourceT $ runConduit $
            filterBin lo hi res (zip (map fst chrSize) rc) .|
            sinkFileBedGzip output
        return $ (fl, emptyFile & location .~ output, n) )
  where
    res = 5000

-- | Make the read count matrix.
mkCellByBinMat :: (Elem 'Gzip tags ~ 'True, SCATACSeqConfig config)
               => SCATACSeq S (File tags 'Bed, File tags 'Bed, Int)
               -> ReaderT config IO (SCATACSeq S (File tags 'Other))
mkCellByBinMat input = do
    dir <- asks ((<> "/ReadCount") . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_readcount.txt.gz" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . (\(tagFl, regionFl, nCell) -> do
        regions <- runResourceT $ runConduit $
            streamBedGzip (regionFl^.location) .| sinkList :: IO [BED3]
        let bedTree = bedToTree undefined $ zip regions [0::Int ..]
            nBin = length regions
            header = B.pack $ printf "Sparse matrix: %d x %d" nCell nBin
        runResourceT $ runConduit $ streamBedGzip (tagFl^.location) .|
            groupBy ((==) `on` (^.name)) .|
            mapC (tagCountPerCell bedTree nBin) .|
            (yield header >> mapC (uncurry showSparseVector)) .|
            unlinesAsciiC .| gzip .| sinkFile output
        return $ emptyFile & location .~ output )

-- | Divide the genome into bins and count the tags.
binCount :: Elem 'Gzip tags ~ 'True
         => File tags 'Bed    -- ^ Tags
         -> [(B.ByteString, Int)]  -- ^ chromosome sizes
         -> Int  -- ^ resolution
         -> IO (Int, [U.Vector Int])   -- ^ Cell number and count for each chromosome
binCount input chrSize res = runResourceT $ runConduit $
    inputWithGroupCell .| zipSinks lengthC (mapC countTags .| fold)
  where
    inputWithGroupCell = streamBedGzip (input^.location) .|
        groupBy ((==) `on` (^.name))
    countTags beds = map (U.map (\x -> if x > 0 then 1 else x)) $
        fst $ runST $ runConduit $ yieldMany beds .| countTagsBinBed res chrs
    chrs = map (\(chr, n) -> BED3 chr 0 n) chrSize
    fold = await >>= \case
        Nothing -> return []
        Just v -> foldlC f v
    f x y = force $ zipWith (U.zipWith (+)) x y
{-# INLINE binCount #-}

filterBin :: Monad m
          => Double
          -> Double
          -> Int
          -> [(B.ByteString, U.Vector Int)] -> ConduitT i BED m ()
filterBin lo hi res rc = forM_ rc $ \(chr, vec) -> flip U.imapM_ vec $ \i x -> 
    when (fromIntegral x > lo && fromIntegral x < hi) $
        yield $ BED chr (i * res) ((i+1) * res) Nothing (Just x) Nothing
{-# INLINE filterBin #-}

tagCountPerCell :: BEDTree Int   -- ^ Regions
                -> Int           -- ^ length of vector
                -> [BED]         -- ^ Tags
                -> (B.ByteString, U.Vector Int)
tagCountPerCell regions n input = (cellBc, rc)
  where
    rc = U.create $ do
        vec <- UM.replicate n 0
        forM_ input $ \x -> forM_ (IM.elems $ intersecting regions $ toSite x) $
            UM.unsafeModify vec (+1)
        return vec
    cellBc = fromJust $ head input^.name
    toSite bed = case bed^.strand of
        Just False -> BED3 (bed^.chrom) (bed^.chromEnd - 1) (bed^.chromEnd)
        _ -> BED3 (bed^.chrom) (bed^.chromStart) (bed^.chromStart + 1)
{-# INLINE tagCountPerCell #-}

showSparseVector :: B.ByteString -> U.Vector Int -> B.ByteString
showSparseVector cellBc = B.intercalate "\t" . (cellBc:) . map f . U.toList .
    U.imapMaybe (\i v -> if v /= 0 then Just (i,v) else Nothing)
  where
    f (i,v) = fromJust (packDecimal i) <> "," <> fromJust (packDecimal v)
{-# INLINE showSparseVector #-}

mergeMatrix :: ( Elem 'Gzip tags1 ~ 'True
               , Elem 'Gzip tags2 ~ 'True
               , SCATACSeqConfig config )
            => [SCATACSeq S (File tags1 'Bed, File tags2 'Other)]
            -> ReaderT config IO (Maybe (File '[Gzip] 'Other))
mergeMatrix [] = return Nothing
mergeMatrix inputs = do
    dir <- asks ((<> "/ReadCount") . _scatacseq_output_dir) >>= getPath
    let output = dir ++ "/Merged_readcount.txt.gz"
    liftIO $ do
        indices <- getIndices $ inputs^..folded.replicates._2.files._1.location
        mats <- forM inputs $ \input -> do
            let e = B.pack $ T.unpack $ input^.eid
                (idxFl, matFl) = input^.replicates._2.files
            idxMap <- fmap (mkIdxMap indices) $ runResourceT $ runConduit $
                streamBedGzip (idxFl^.location) .|
                mapC (\x -> ((x::BED3)^.chrom, x^.chromStart)) .| sinkList
            mat <- mkSpMatrix readInt $ matFl^.location
            return ( e, mat
                { _decoder = \x -> _decoder mat x .| mapC (changeIdx idxMap)
                , _num_col = M.size indices } )
        merge output mats
    return $ Just $ location .~ output $ emptyFile
  where
    merge :: FilePath -> [(B.ByteString, SpMatrix Int)] -> IO ()
    merge output ms
        | any (/=nBin) (map (_num_col . snd) ms) = error "Column unmatched!"
        | otherwise = runResourceT $ runConduit $ source .|
            (yield header >> mapC (encodeRowWith (fromJust . packDecimal))) .|
            unlinesAsciiC .| gzip .| sinkFile output
      where
        source = forM_ ms $ \(sample, m) -> streamRows m .| 
            mapC (first (\x -> sample <> "_" <> x))
        nCell = foldl' (+) 0 $ map (_num_row . snd) ms
        nBin = _num_col $ snd $ head ms
        header = B.pack $ printf "Sparse matrix: %d x %d" nCell nBin

changeIdx :: U.Vector Int -> Row a -> Row a
changeIdx idxMap = second $ sortBy (comparing fst) . map (first (idxMap U.!))
{-# INLINE changeIdx #-}

-- | A map from old index to new index.
mkIdxMap :: M.HashMap (B.ByteString, Int) Int
          -> [(B.ByteString, Int)]
          -> U.Vector Int
mkIdxMap newIdx oldIdx = U.fromList $ map f oldIdx
  where
    f k = M.lookupDefault undefined k newIdx
{-# INLINE mkIdxMap #-}

getIndices :: [FilePath]   -- ^ A list of files containg the indices
           -> IO (M.HashMap (B.ByteString, Int) Int)
getIndices idxFls = do
    idx <- runResourceT $ runConduit $
        mapM_ streamBedGzip idxFls .| foldlC f S.empty
    return $ M.fromList $ zip (S.toList idx) [0..]
  where
    f :: S.HashSet (B.ByteString, Int) -> BED3 -> S.HashSet (B.ByteString, Int) 
    f m b = S.insert (b^.chrom, b^.chromStart) m
{-# INLINE getIndices #-}


--------------------------------------------------------------------------------
-- Cell by Gene Matrix
--------------------------------------------------------------------------------

estimateExpr :: SCATACSeqConfig config
             => (SCATACSeq S (B.ByteString, File '[Gzip] 'Bed))
             -> ReaderT config IO
                (SCATACSeq S (B.ByteString, File '[GeneQuant] 'Tsv))
estimateExpr input = do
    genes <- asks _scatacseq_annotation >>= liftIO . readGenes . fromJust
    let idRep = asDir $ "/GeneQuant/" <> T.unpack (input^.eid) <>
            "_rep" <> show (input^.replicates._1)
    dir <- asks _scatacseq_output_dir >>= getPath . (<> idRep)
    input & replicates.traverse.files %%~ liftIO . ( \(nm, fl) -> do
        let output = dir ++ "/" ++ B.unpack nm ++ ".tsv" 
        counts <- runResourceT $ runConduit $
            streamBedGzip (fl^.location) .| mkGeneCount genes
        B.writeFile output $ B.unlines $
            map (\(n, c) -> n <> "\t" <> toShortest c) counts
        return (nm, location .~ output $ emptyFile) )

-- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6385419/
mkGeneCount :: PrimMonad m
            => [Gene]
            -> ConduitT BED o m [(B.ByteString, Double)]
mkGeneCount genes = zip (map (original . geneName) genes) . U.toList <$>
    rpkmBed (map getTss genes)
  where
    getTss Gene{..} | geneStrand = BED3 geneChrom (geneLeft - 1000) (geneLeft + 1000)
                    | otherwise = BED3 geneChrom (geneRight - 1000) (geneRight + 1000)
{-# INLINE mkGeneCount #-}

-- | Combine expression data into a table and output
mkExprTable :: SCATACSeqConfig config
            => [SCATACSeq S (B.ByteString, File '[GeneQuant] 'Tsv)]
            -> ReaderT config IO
                [SCATACSeq S (File '[GeneQuant] 'Tsv)]
mkExprTable = mapM mkTable . concatMap split . mergeExp
  where
    mkTable input = do
        results <- liftIO $ readExpr input
        let idRep = asDir $ "/GeneQuant/" <> T.unpack (input^.eid) <>
                "_rep" <> show (input^.replicates._1)
        dir <- asks _scatacseq_output_dir >>= getPath . (<> idRep)
        let output = dir ++ "/" ++ "expression_profile.tsv"
            (expNames, values) = unzip $ M.toList $
                fmap (map (second average) . combine) $ M.fromListWith (++) $ results
        liftIO $ B.writeFile output $ B.unlines $
            B.intercalate "\t" ("Name" : expNames) :
            map (\(x,xs) -> B.intercalate "\t" $ original x : map toShortest xs)
                (combine values)
        return $ input & replicates._2.files .~ (location .~ output $ emptyFile)

--------------------------------------------------------------------------------
-- Auxiliary functions
--------------------------------------------------------------------------------

readExpr :: SCATACSeq S [(B.ByteString, File '[GeneQuant] 'Tsv)]
         -> IO [(B.ByteString, [[(CI B.ByteString, Double)]])]
readExpr input = forM (input^.replicates._2.files) $ \(nm, fl) -> do
    c <- B.readFile $ fl^.location
    let result = map (\xs ->
            let fs = B.split '\t' xs in (mk $ head fs, readDouble $ fs!!1)) $
            tail $ B.lines c
    return (nm, [result])
{-# INLINE readExpr #-}

combine :: [[(CI B.ByteString, Double)]] -> [(CI B.ByteString, [Double])]
combine xs = flip map names $ \nm -> (nm, map (M.lookupDefault 0.01 nm) xs')
  where
    names = nubSort $ concatMap (fst . unzip) xs
    xs' = map (fmap average . M.fromListWith (++) . map (second return)) xs
{-# INLINE combine #-}

average :: [Double] -> Double
average [a,b]   = (a + b) / 2
average [a,b,c] = (a + b + c) / 3
average xs      = foldl1' (+) xs / fromIntegral (length xs)
{-# INLINE average #-}