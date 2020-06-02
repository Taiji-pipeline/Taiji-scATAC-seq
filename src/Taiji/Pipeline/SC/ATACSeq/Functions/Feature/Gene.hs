{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE GADTs #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Feature.Gene
    ( mkExprTable
    , writePromoters
    , mkCellByGene
    , GeneAccDef(..)
    ) where

import qualified Data.ByteString.Char8 as B
import qualified Data.HashMap.Strict                  as M
import Control.Arrow (second)
import Data.Char (toUpper)
import Bio.Data.Bed.Types
import Bio.Data.Bed
import Data.Singletons.Prelude (Elem)
import qualified Data.Matrix            as Mat
import qualified Data.Text as T
import Bio.RealWorld.GENCODE (Gene(..), Transcript(..), TranscriptType(..))
import           Data.CaseInsensitive  (original)
import           Bio.Pipeline.Utils
import qualified Data.Vector as V
import Data.Conduit.Zlib (gzip)
import qualified Data.IntervalMap.Strict as IM

import Taiji.Prelude
import Taiji.Utils
import Taiji.Pipeline.SC.ATACSeq.Functions.Utils
import Taiji.Pipeline.SC.ATACSeq.Types
import qualified Taiji.Utils.DataFrame as DF

mkCellByGene :: (Elem 'Gzip tags ~ 'True, SCATACSeqConfig config)
             => FilePath
             -> (SCATACSeq S (File tags 'Bed, Int), File '[] 'Tsv)
             -> ReaderT config IO (SCATACSeq S (File '[Gzip] 'Other))
mkCellByGene prefix (input, promoters) = do
    dir <- asks ((<> asDir prefix) . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_gene.mat.gz" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . ( \(fl, nCell) -> do
        regions <- map (map readBED3 . B.split ',' . last . B.split '\t') . B.lines <$>
            B.readFile (promoters^.location)
        runResourceT $ runConduit $ streamBedGzip (fl^.location) .|
            groupCells .| mkMat nCell regions .| sinkFile output
        return $ emptyFile & location .~ output )
  where
    mkMat :: (PrimMonad m, MonadThrow m)
          => Int   -- ^ the number of cells
          -> [[BED3]]    -- ^ a list of regions
          -> ConduitT (B.ByteString, [BED]) B.ByteString m ()
    mkMat nCell regions = source .| unlinesAsciiC .| gzip
      where
        source = yield header >> mapC
            (encodeRowWith (fromJust . packDecimal) . second (countEachCell bedTree))
        bedTree = bedToTree (++) $ concat $
            zipWith (\i xs -> zipWith (\j x -> (x, [(i,j)])) [0..] xs) [0..] regions
        header = B.pack $ printf "Sparse matrix: %d x %d" nCell (length regions)
        countEachCell :: BEDTree [(Int, Int)] -> [BED] -> [(Int, Int)]
        countEachCell beds = M.toList . M.fromListWith max .
            map (\((i,_), x) -> (i,x)) . M.toList . foldl' f M.empty
          where
            f m bed = foldl' (\x k -> M.insertWith (+) k (1::Int) x) m $
                concat $ intersecting beds query
              where
                query = case bed^.strand of
                    Just False -> BED3 (bed^.chrom) (bed^.chromEnd - 1) (bed^.chromEnd)
                    _ -> BED3 (bed^.chrom) (bed^.chromStart) (bed^.chromStart + 1)

data GeneAccDef = PromoterOnly -- ^ -/+ 1000 around TSS
                | PromoterPlusGeneBody   -- ^ Gene body plus upstream 2000

writePromoters :: SCATACSeqConfig config 
               => GeneAccDef 
               -> ReaderT config IO (File '[] 'Tsv)
writePromoters def = do
    dir <- asks ((<> "/Feature/Gene/") . _scatacseq_output_dir) >>= getPath
    let output = dir <> "promoters.tsv"
    genes <- asks _scatacseq_annotation >>= liftIO . readGenesValidated . fromJust
    let promoters = concatMap getPromoter genes
    liftIO $ do
        let f xs = fromJust (head xs^.name) <> "\t" <>
                B.intercalate "," (map showBED xs)
        B.writeFile output $ B.unlines $ map f $ groupBy ((==) `on` (^.name)) $
            sortBy (comparing (^.name)) promoters
    return $ location .~ output $ emptyFile
  where
    getPromoter Gene{..} = case def of
        PromoterOnly -> map f $ keepCoding geneTranscripts
        PromoterPlusGeneBody -> case geneStrand of
            True -> [BED geneChrom (max 0 $ geneLeft - 2000)
                geneRight nm Nothing (Just geneStrand)]
            False -> [BED geneChrom geneLeft
                (geneRight + 2000) nm Nothing (Just geneStrand)]
      where
        keepCoding xs =
            let xs' = filter ((==Coding) . transType) xs
            in if null xs' then xs else xs'
        nm = Just $ B.map toUpper $ original geneName
        f Transcript{..}
            | transStrand = BED geneChrom (max 0 $ transLeft - 1000)
                (transLeft + 1000) nm Nothing (Just transStrand)
            | otherwise = BED geneChrom (max 0 $ transRight - 1000)
                (transRight + 1000) nm Nothing (Just transStrand)

-- | Combine expression data into a table and output
mkExprTable :: SCATACSeqConfig config
            => FilePath
            -> ( Maybe (File '[] 'Tsv) -- ^ genes
               , Maybe (FilePath, FilePath, FilePath) )
            -> ReaderT config IO (Maybe (File '[GeneQuant] 'Tsv))
mkExprTable prefix (Just fl, Just (accFl,_,_)) = do
    dir <- asks _scatacseq_output_dir >>= getPath . (<> asDir prefix)
    liftIO $ do
        genes <- map
            ((\[nm, x] -> (nm, map readBED3 $ B.split ',' x)) . B.split '\t') .
            B.lines <$> B.readFile (fl^.location)
        let output = dir ++ "/gene_accessibility.tsv"
        acc <- DF.readTable accFl
        DF.writeTable output (T.pack . show) $ mkTable genes acc
        return $ Just $ location .~ output $ emptyFile
mkExprTable _ _ = return Nothing

mkTable :: [(B.ByteString, [BED3])]  -- ^ Promoters
        -> DF.DataFrame Double   -- ^ accessibility
        -> DF.DataFrame Double
mkTable promoters acc = DF.fromMatrix rows (DF.colNames acc) $ Mat.fromRows vals
  where
    (rows, vals) = unzip $ mapMaybe f promoters
    f (nm, pro) = do
        vec <- foldl1'' (V.zipWith max) $ flip mapMaybe pro $
            foldl1'' (V.zipWith (+)) . map (acc `DF.rindex`) . IM.elems .
                intersecting peaks
        return (T.pack $ B.unpack nm, vec)
    peaks = bedToTree undefined $ zip
        (map (readBED3 . B.pack . T.unpack) $ DF.rowNames acc) [0 :: Int ..]
    foldl1'' _ [] = Nothing
    foldl1'' g xs = Just $ foldl1' g xs
{-# INLINE mkTable #-}

showBED :: BED -> B.ByteString
showBED x = (x^.chrom) <> ":" <>
    (fromJust $ packDecimal $ x^.chromStart) <> "-" <>
    (fromJust $ packDecimal $ x^.chromEnd)
{-# INLINE showBED #-}

readBED3 :: B.ByteString -> BED3
readBED3 x = let [chr, x'] = B.split ':' x
                 [s, e] = B.split '-' x'
             in asBed chr (readInt s) $ readInt e
{-# INLINE readBED3 #-}