{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE GADTs #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Feature
    ( mkCutSiteIndex
    , getBins
    , mkWindowMat

    , mkPeakMat
    , findPeaks
    , mergePeaks

    , mergeFeatMatrix
    , estimateExpr
    , mkExprTable
    ) where

import qualified Data.ByteString.Char8 as B
import Data.Conduit.List (groupBy)
import qualified Data.HashMap.Strict                  as M
import Bio.Data.Bed.Types
import Bio.Data.Bed
import Bio.RealWorld.GENCODE (readGenes, Gene(..))
import Data.Double.Conversion.ByteString (toShortest)
import Bio.Utils.Misc (readDouble)
import           Data.CaseInsensitive  (mk, original, CI)
import qualified Data.Vector.Unboxed as U
import Data.Singletons.Prelude (Elem)
import           Bio.Pipeline.Utils
import Control.Arrow (second)
import           Data.List.Ordered                    (nubSort)
import qualified Data.Text as T
import Bio.Data.Bed.Utils
import Bio.Seq.IO (withGenome, getChrSizes)

import Taiji.Prelude hiding (groupBy)
import Taiji.Pipeline.SC.ATACSeq.Types
import Taiji.Pipeline.SC.ATACSeq.Functions.Utils
import Taiji.Pipeline.SC.ATACSeq.Functions.Feature.Window
import Taiji.Pipeline.SC.ATACSeq.Functions.Feature.Peak

mergeFeatMatrix :: ( Elem 'Gzip tags1 ~ 'True
                   , Elem 'Gzip tags2 ~ 'True
                   , SCATACSeqConfig config )
                => FilePath
                -> [SCATACSeq S (File tags1 'Bed, File tags2 'Other)]
                -> ReaderT config IO [SCATACSeq S (File '[Gzip] 'Other)]
mergeFeatMatrix _ [] = return []
mergeFeatMatrix filename inputs = do
    dir <- asks ((<> "/Merged") . _scatacseq_output_dir) >>= getPath
    let output = dir ++ "/" ++ filename
    liftIO $ runResourceT $ runConduit $ mergeMatrix inputs' .| sinkFile output
    return $ return $ (head inputs & eid .~ "Merged") & replicates._2.files .~
        (location .~ output $ emptyFile)
  where
    inputs' = map (\x -> (B.pack $ T.unpack $ x^.eid, x^.replicates._2.files)) inputs

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