{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE DataKinds #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Feature.Gene
    ( estimateExpr
    , mkExprTable
    ) where

import qualified Data.ByteString.Char8 as B
import qualified Data.HashMap.Strict                  as M
import Bio.Data.Bed.Types
import Bio.Data.Bed
import Bio.Data.Bed.Utils (rpkmBed)
import Bio.RealWorld.GENCODE (readGenes, Gene(..))
import Data.Double.Conversion.ByteString (toShortest)
import Bio.Utils.Misc (readDouble)
import           Data.CaseInsensitive  (mk, original, CI)
import qualified Data.Vector.Unboxed as U
import           Bio.Pipeline.Utils
import Control.Arrow (second)
import           Data.List.Ordered                    (nubSort)

import Taiji.Prelude hiding (groupBy)
import Taiji.Pipeline.SC.ATACSeq.Types


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

{-
getTSS :: FilePath -> IO (BEDTree (Int, Bool))
getTSS fl = do
    genes <- readGenes fl
    map (original . geneName)
    
    fmap (bedToTree const . concatMap fn) . readGenes
  where
    fn Gene{..} = map g $ nubSort tss
      where
        g x = (BED3 geneChrom (x - 1000) (x + 1000), (x, geneStrand))
        tss | geneStrand = geneLeft : map fst geneTranscripts
            | otherwise = geneRight : map snd geneTranscripts
-}
 

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