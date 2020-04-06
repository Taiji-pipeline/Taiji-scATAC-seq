{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE GADTs #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Feature.Gene
    ( mkExprTable
    , writePromoters
    ) where

import qualified Data.ByteString.Char8 as B
import qualified Data.HashMap.Strict                  as M
import Control.Arrow (first, second)
import Data.Char (toUpper)
import Bio.Data.Bed.Types
import Bio.Data.Bed
import Data.List.Ordered (nubSort)
import Data.Singletons.Prelude (Elem)
import qualified Data.Matrix            as Mat
import qualified Data.Text as T
import Bio.RealWorld.GENCODE (readGenes, Gene(..), Transcript(..), TranscriptType(..))
import           Data.CaseInsensitive  (original)
import           Bio.Pipeline.Utils
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector as V
import Data.Conduit.Zlib (gzip)
import qualified Data.IntervalMap.Strict as IM

import Taiji.Prelude
import Taiji.Utils
import Taiji.Pipeline.SC.ATACSeq.Functions.Utils
import Taiji.Pipeline.SC.ATACSeq.Types
import qualified Taiji.Utils.DataFrame as DF

{-
mkCellByGene :: (Elem 'Gzip tags ~ 'True, SCATACSeqConfig config)
             => FilePath
             -> (SCATACSeq S (File tags 'Bed, Int), (File '[Gzip] 'Bed, File '[] 'Tsv))
             -> ReaderT config IO (SCATACSeq S (File '[Gzip] 'Other))
mkCellByGene prefix (input, (promoters, geneFl)) = do
    dir <- asks ((<> asDir prefix) . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_gene.mat.gz" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . ( \(fl, nCell) -> do
        regions <- runResourceT $ runConduit $
            streamBedGzip (promoters^.location) .| sinkList
        genes <- fmap B.lines $ B.readFile $ geneFl^.location
        let geneIdx = M.fromList $ zip genes [0..]
        runResourceT $ runConduit $ streamBedGzip (fl^.location) .|
            groupCells .| mkMat nCell regions geneIdx .| sinkFile output
        return $ emptyFile & location .~ output )
  where
    mkMat :: (PrimMonad m, MonadThrow m)
          => Int   -- ^ the number of cells
          -> [BED]    -- ^ a list of regions
          -> M.HashMap B.ByteString Int
          -> ConduitT (B.ByteString, [BED]) B.ByteString m ()
    mkMat nCell regions geneIdx = source .| unlinesAsciiC .| gzip
      where
        names = V.fromList $ map (fromJust . (^.name)) regions
        source = yield header >> mapC
            (encodeRowWith (fromJust . packDecimal) . second (countEachCell bedTree))
        bedTree = bedToTree (++) $ zip regions $ map return [0::Int ..]
        nBin = M.size geneIdx
        header = B.pack $ printf "Sparse matrix: %d x %d" nCell nBin
        countEachCell :: BEDTree [Int] -> [BED] -> [(Int, Int)]
        countEachCell beds = M.toList . M.fromListWith max .
            map (first findIdx) . M.toList . foldl' f M.empty
          where
            findIdx i = M.lookupDefault undefined (names V.! i) geneIdx
            f m bed = foldl' (\x k -> M.insertWith (+) k (1::Int) x) m $
                concat $ intersecting beds query
              where
                query = case bed^.strand of
                    Just False -> BED3 (bed^.chrom) (bed^.chromEnd - 1) (bed^.chromEnd)
                    _ -> BED3 (bed^.chrom) (bed^.chromStart) (bed^.chromStart + 1)
-}

writePromoters :: SCATACSeqConfig config 
               => ReaderT config IO (File '[] 'Tsv)
writePromoters = do
    dir <- asks ((<> "/Feature/Gene/") . _scatacseq_output_dir) >>= getPath
    let output = dir <> "promoters.tsv"
    genes <- asks _scatacseq_annotation >>= liftIO . readGenes . fromJust
    let promoters = concatMap getPromoter genes
    liftIO $ do
        let f xs = fromJust (head xs^.name) <> "\t" <>
                B.intercalate "," (map showBED xs)
        B.writeFile output $ B.unlines $ map f $ groupBy ((==) `on` (^.name)) $ sortBy (comparing (^.name)) promoters
    return $ location .~ output $ emptyFile
  where
    getPromoter gene = map f $ keepCoding $ geneTranscripts gene
      where
        keepCoding xs =
            let xs' = filter ((==Coding) . transType) xs
            in if null xs' then xs else xs'
        nm = Just $ B.map toUpper $ original $ geneName gene
        f Transcript{..}
            | transStrand = BED (geneChrom gene) (max 0 $ transLeft - 1000)
                (transLeft + 1000) nm Nothing (Just transStrand)
            | otherwise = BED (geneChrom gene) (max 0 $ transRight - 1000)
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