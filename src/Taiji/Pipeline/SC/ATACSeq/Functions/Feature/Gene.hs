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
import Bio.Data.Bed.Types
import Bio.Data.Bed
import Data.Singletons.Prelude (Elem)
import qualified Data.Text as T
import Bio.RealWorld.GENCODE (readGenes, Gene(..), Transcript(..))
import           Data.CaseInsensitive  (original)
import           Bio.Pipeline.Utils
import qualified Data.Vector.Unboxed as U

import Taiji.Prelude
import Taiji.Pipeline.SC.ATACSeq.Functions.Utils
import Taiji.Pipeline.SC.ATACSeq.Types
import qualified Taiji.Utils.DataFrame as DF

mkCellByGene :: (Elem 'Gzip tags ~ 'True, SCATACSeqConfig config)
             => FilePath
             -> (SCATACSeq S (File tags 'Bed, Int), File '[Gzip] 'Bed)
             -> ReaderT config IO (SCATACSeq S (File '[Gzip] 'Other))
mkCellByGene prefix (input, promoters) = do
    dir <- asks ((<> asDir prefix) . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_cell_by_transcript.mat.gz" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . ( \(fl, nCell) -> do
        regions <- fmap (map return) $ runResourceT $ runConduit $
            streamBedGzip (promoters^.location) .| sinkList
        runResourceT $ runConduit $ streamBedGzip (fl^.location) .|
            groupCells .| mkFeatMat nCell regions .| sinkFile output
        return $ emptyFile & location .~ output )

writePromoters :: SCATACSeqConfig config 
               => ReaderT config IO (File '[Gzip] 'Bed)
writePromoters = do
    dir <- asks ((<> "/Feature/Gene/") . _scatacseq_output_dir) >>= getPath
    let output = dir <> "promoters.bed.gz"
    genes <- asks _scatacseq_annotation >>= liftIO . readGenes . fromJust
    liftIO $ runResourceT $ runConduit $
        yieldMany (concatMap getPromoter genes) .| sinkFileBedGzip output
    return $ location .~ output $ emptyFile
  where
    getPromoter gene = map f $ geneTranscripts gene
      where
        f Transcript{..}
            | transStrand = BED (geneChrom gene) (max 0 $ transLeft - 1000)
                (transLeft + 1000) (Just transId) Nothing (Just transStrand)
            | otherwise = BED (geneChrom gene) (max 0 $ transRight - 1000)
                (transRight + 1000) (Just transId) Nothing (Just transStrand)

-- | Combine expression data into a table and output
mkExprTable :: SCATACSeqConfig config
            => FilePath
            -> ( File '[Gzip] 'Bed -- ^ genes
               , [SCATACSeq S (File '[Gzip] 'Other)] )
            -> ReaderT config IO (Maybe (File '[GeneQuant] 'Tsv))
mkExprTable prefix (fl, inputs) = do
    dir <- asks _scatacseq_output_dir >>= getPath . (<> asDir prefix)
    idToGene <- asks _scatacseq_annotation >>= liftIO . readGenes . fromJust >>=
        return . M.fromList . concatMap f
    liftIO $ do
        promoters <- runResourceT $ runConduit $ streamBedGzip (fl^.location) .| sinkList :: IO [BED]
        let geneNames = map (\x -> M.lookupDefault undefined (fromJust $ x^.name) idToGene) promoters
            output = dir ++ "/gene_accessibility.tsv"

        mat <- fmap transpose $ forM inputs $ \input -> fmap U.toList $
            computeRAS $ input^.replicates._2.files.location 
        let (genes, vals) = unzip $ map combine $ groupBy ((==) `on` fst) $
                sortBy (comparing fst) $ zip geneNames mat
        DF.writeTable output (T.pack . show) $
            DF.mkDataFrame (map (T.pack . B.unpack) genes) (map (^.eid) inputs) vals
        return $ Just $ location .~ output $ emptyFile
  where
    f Gene{..} = zip (map transId geneTranscripts) $
        repeat $ original geneName
    combine xs = (head gene, foldl1' (zipWith max) vals)
      where
        (gene, vals) = unzip xs