{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE GADTs #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Feature.Gene
    ( mkExprTable
    , writePromoters
    , mkCellByGene
    ) where

import qualified Data.ByteString.Char8 as B
import qualified Data.HashMap.Strict                  as M
import Control.Arrow (first, second)
import Data.Char (toUpper)
import Bio.Data.Bed.Types
import Bio.Data.Bed
import Data.List.Ordered (nubSort)
import Data.Singletons.Prelude (Elem)
import qualified Data.Text as T
import Bio.RealWorld.GENCODE (readGenes, Gene(..), Transcript(..))
import           Data.CaseInsensitive  (original)
import           Bio.Pipeline.Utils
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector as V
import Data.Conduit.Zlib (gzip)

import Taiji.Prelude
import Taiji.Utils
import Taiji.Pipeline.SC.ATACSeq.Functions.Utils
import Taiji.Pipeline.SC.ATACSeq.Types
import qualified Taiji.Utils.DataFrame as DF

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

writePromoters :: SCATACSeqConfig config 
               => ReaderT config IO (File '[Gzip] 'Bed, File '[] 'Tsv)
writePromoters = do
    dir <- asks ((<> "/Feature/Gene/") . _scatacseq_output_dir) >>= getPath
    let output = dir <> "promoters.bed.gz"
        geneOutput = dir <> "gene_idx.txt"
    genes <- asks _scatacseq_annotation >>= liftIO . readGenes . fromJust
    let promoters = concatMap getPromoter genes
    liftIO $ do
        runResourceT $ runConduit $
            yieldMany promoters .| sinkFileBedGzip output
        B.writeFile geneOutput $ B.unlines $ nubSort $ map (fromJust . (^.name)) promoters
    return ( location .~ output $ emptyFile
           , location .~ geneOutput $ emptyFile )
  where
    getPromoter gene = map f $ geneTranscripts gene
      where
        nm = Just $ B.map toUpper $ original $ geneName gene
        f Transcript{..}
            | transStrand = BED (geneChrom gene) (max 0 $ transLeft - 1000)
                (transLeft + 1000) nm Nothing (Just transStrand)
            | otherwise = BED (geneChrom gene) (max 0 $ transRight - 1000)
                (transRight + 1000) nm Nothing (Just transStrand)

-- | Combine expression data into a table and output
mkExprTable :: SCATACSeqConfig config
            => FilePath
            -> ( Maybe (a, File '[] 'Tsv) -- ^ genes
               , [SCATACSeq S (File '[Gzip] 'Other)] )
            -> ReaderT config IO (Maybe (File '[GeneQuant] 'Tsv))
mkExprTable prefix (Just (_, fl), inputs) = do
    dir <- asks _scatacseq_output_dir >>= getPath . (<> asDir prefix)
    liftIO $ do
        geneNames <- fmap B.lines $ B.readFile $ fl^.location
        let output = dir ++ "/gene_accessibility.tsv"
            output2 = dir ++ "/gene_specificity.tsv"
        mat <- fmap transpose $ forM inputs $ \input -> fmap U.toList $
            computeRAS $ input^.replicates._2.files.location 
        let (genes, vals) = unzip $ map combine $ groupBy ((==) `on` fst) $
                sortBy (comparing fst) $ zip geneNames mat
            ras = DF.mkDataFrame (map (T.pack . B.unpack) genes)
                (map (^.eid) inputs) vals
            ss = computeSS ras
        DF.writeTable output (T.pack . show) ras
        DF.writeTable output2 (T.pack . show) ss
        return $ Just $ location .~ output $ emptyFile
  where
    combine xs = (head gene, foldl1' (zipWith max) vals)
      where
        (gene, vals) = unzip xs
mkExprTable _ _ = return Nothing