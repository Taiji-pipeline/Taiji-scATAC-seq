{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE GADTs #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Feature.Gene
    ( estimateExpr
    , mkExprTable
    , getGeneNames
    , mkCellByGene 
    , mergeCellByGeneMatrix
    ) where

import qualified Data.ByteString.Char8 as B
import qualified Data.HashMap.Strict                  as M
import Data.ByteString.Lex.Integral (packDecimal)
import Bio.Data.Bed.Types
import Bio.Data.Bed
import Data.Singletons.Prelude (Elem)
import qualified Data.Text as T
import Control.Arrow (first)
import Bio.Data.Bed.Utils (rpkmBed)
import Data.Conduit.Zlib (gzip)
import Bio.RealWorld.GENCODE (readGenes, Gene(..))
import Data.Double.Conversion.ByteString (toShortest)
import Bio.Utils.Misc (readDouble, readInt)
import           Data.CaseInsensitive  (mk, original, CI)
import qualified Data.Vector.Unboxed as U
import           Bio.Pipeline.Utils
import Control.Arrow (second)
import           Data.List.Ordered                    (nubSort)

import Taiji.Prelude hiding (groupBy)
import Taiji.Pipeline.SC.ATACSeq.Functions.Utils
import Taiji.Pipeline.SC.ATACSeq.Types

mkCellByGene :: (Elem 'Gzip tags ~ 'True, SCATACSeqConfig config)
             => FilePath
             -> (SCATACSeq S (File tags 'Bed, a, Int), FilePath)
             -> ReaderT config IO (SCATACSeq S (FilePath, File '[Gzip] 'Other))
mkCellByGene prefix (input, genes) = do
    dir <- asks ((<> asDir prefix) . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_cell_by_gene.mat.gz" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . ( \(fl,_,nCell) -> do
        regions <- map snd <$> readTSS genes
        runResourceT $ runConduit $ streamBedGzip (fl^.location) .|
            groupCells .| mkFeatMat nCell regions .| sinkFile output
        return $ (genes, emptyFile & location .~ output) )

mergeCellByGeneMatrix :: SCATACSeqConfig config
                      => [SCATACSeq S (a, File '[Gzip] 'Other)]
                      -> ReaderT config IO [SCATACSeq S (a, File '[Gzip] 'Other)]
mergeCellByGeneMatrix inputs = do
    dir <- asks ((<> "/Feature/Gene/") . _scatacseq_output_dir) >>= getPath
    let output = dir <> "/Merged_cell_by_gene.mat.gz"
    mats <- liftIO $ forM inputs $ \input -> do
        let nm = B.pack $ T.unpack $ input^.eid
        mat <- mkSpMatrix id $ input^.replicates._2.files._2.location
        return (nm, mat)
    liftIO $ runResourceT $ runConduit $ merge mats .| sinkFile output
    return $ return $ head inputs & eid .~ "Merged" &
        replicates._2.files._2.location .~ output
  where
    merge ms
        | any (/=nBin) (map (_num_col . snd) ms) = error "Column unmatched!"
        | otherwise = source .| (yield header >> mapC (encodeRowWith id)) .|
            unlinesAsciiC .| gzip
      where
        source = forM_ ms $ \(sample, m) -> streamRows m .| 
            mapC (first (\x -> sample <> "+" <> x))
        nCell = foldl' (+) 0 $ map (_num_row . snd) ms
        nBin = _num_col $ snd $ head ms
        header = B.pack $ printf "Sparse matrix: %d x %d" nCell nBin

getGeneNames :: SCATACSeqConfig config
             => ReaderT config IO FilePath
getGeneNames = do
    dir <- asks ((<> "/Feature/Gene/") . _scatacseq_output_dir) >>= getPath
    let output = dir <> "gene_name_idx.tsv"
    tss <- asks _scatacseq_annotation >>= liftIO . getTSS . fromJust
    liftIO $ writeTSS output tss
    return output

writeTSS :: FilePath -> [(B.ByteString, [BED3])] -> IO ()
writeTSS output xs = B.writeFile output $ B.unlines $ flip map xs $
    \(nm, regions) -> nm <> "\t" <> B.intercalate "," (map f regions)
  where
    f b = b^.chrom <> ":" <> (fromJust $ packDecimal $ b^.chromStart) <> "-" <>
        (fromJust $ packDecimal $ b^.chromEnd) 

readTSS :: FilePath -> IO [(B.ByteString, [BED3])]
readTSS input = map f . B.lines <$> B.readFile input
  where
    f x = let [nm, regions] = B.split '\t' x
          in (nm, map g $ B.split ',' regions)
    g x = let [chr, y] = B.split ':' x
              [s,e] = B.split '-' y
          in asBed chr (readInt s) $ readInt e

getTSS :: FilePath -> IO [(B.ByteString, [BED3])]
getTSS fl = do
    genes <- nubSort . concatMap fn <$> readGenes fl
    return $ flip map (M.toList $ M.fromListWith (++) genes) $ \(gene, regions) ->
        (original gene, runIdentity $ runConduit $ mergeBed regions .| sinkList)
  where
    fn Gene{..} = map g $ nubSort tss
      where
        g x = (geneName, [BED3 geneChrom (max 0 $ x - 1000) (x + 1000)])
        tss | geneStrand = geneLeft : map fst geneTranscripts
            | otherwise = geneRight : map snd geneTranscripts

-- | Estimate the gene expression level using atac-seq counts.
estimateExpr :: SCATACSeqConfig config
             => (B.ByteString, File '[Gzip] 'Bed)
             -> ReaderT config IO (B.ByteString, File '[GeneQuant] 'Tsv)
estimateExpr (nm, fl) = do
    genes <- asks _scatacseq_annotation >>= liftIO . readGenes . fromJust
    dir <- asks _scatacseq_output_dir >>= getPath . (<> "/Feature/Gene/")
    let output = dir ++ "/" ++ B.unpack nm ++ ".tsv" 
    liftIO $ do
        counts <- runResourceT $ runConduit $
            streamBedGzip (fl^.location) .| mkGeneCount genes
        B.writeFile output $ B.unlines $
            map (\(n, c) -> n <> "\t" <> toShortest c) counts
        return (nm, location .~ output $ emptyFile)

-- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6385419/
-- | Count the tags in promoter regions (RPKM).
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
            => FilePath
            -> [(B.ByteString, File '[GeneQuant] 'Tsv)]
            -> ReaderT config IO (Maybe (File '[GeneQuant] 'Tsv))
mkExprTable _ [] = return Nothing
mkExprTable nm input = do
    results <- liftIO $ readExpr input
    dir <- asks _scatacseq_output_dir >>= getPath . (<> "/Feature/Gene/")
    let output = dir ++ nm
        (expNames, values) = unzip $ M.toList $
            fmap (map (second average) . combine) $ M.fromListWith (++) $ results
    liftIO $ B.writeFile output $ B.unlines $
        B.intercalate "\t" ("Name" : expNames) :
        map (\(x,xs) -> B.intercalate "\t" $ original x : map toShortest xs)
            (combine values)
    return $ Just $ location .~ output $ emptyFile

--------------------------------------------------------------------------------
-- Auxiliary functions
--------------------------------------------------------------------------------

readExpr :: [(B.ByteString, File '[GeneQuant] 'Tsv)]
         -> IO [(B.ByteString, [[(CI B.ByteString, Double)]])]
readExpr input = forM input $ \(nm, fl) -> do
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