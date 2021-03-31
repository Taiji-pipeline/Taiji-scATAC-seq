{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE GADTs #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Feature.Gene
    ( mkExprTable
    , combineExprTable
    , writePromoters
    , mkCellByGene
    , GeneAccDef(..)
    ) where

import qualified Data.ByteString.Char8 as B
import qualified Data.HashMap.Strict                  as M
import Data.Binary (decodeFile, encodeFile)
import Control.Arrow (second)
import Data.Char (toUpper)
import Bio.Data.Bed.Types
import Bio.Data.Bed
import Data.Singletons.Prelude (Elem)
import qualified Data.Text as T
import Bio.RealWorld.GENCODE (Gene(..), Transcript(..), TranscriptType(..))
import           Data.CaseInsensitive  (original)
import           Bio.Pipeline.Utils
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Unboxed.Mutable as UM
import Data.Conduit.Zlib (gzip)
import qualified Data.IntervalMap.Strict as IM
import Control.DeepSeq (force)

import Taiji.Prelude
import Taiji.Utils
import Taiji.Pipeline.SC.ATACSeq.Functions.Utils
import Taiji.Pipeline.SC.ATACSeq.Types
import qualified Taiji.Utils.DataFrame as DF

mkCellByGene :: (Elem 'Gzip tags ~ 'True, SCATACSeqConfig config)
             => FilePath
             -> (SCATACSeq S (File tags 'Bed, Int), File '[] 'Bed)
             -> ReaderT config IO (SCATACSeq S (File '[Gzip] 'Other))
mkCellByGene prefix (input, promoters) = do
    dir <- asks ((<> asDir prefix) . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_gene.mat.gz" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . ( \(fl, nCell) -> do
        regions <- readBed $ promoters^.location :: IO [BED3]
        runResourceT $ runConduit $ streamBedGzip (fl^.location) .|
            groupCells .| mkMat nCell regions .| sinkFile output
        return $ emptyFile & location .~ output )
  where
    mkMat :: (PrimMonad m, MonadThrow m)
          => Int   -- ^ the number of cells
          -> [BED3]    -- ^ a list of regions
          -> ConduitT (B.ByteString, [BED]) B.ByteString m ()
    mkMat nCell regions = source .| unlinesAsciiC .| gzip
      where
        source = yield header >> mapC
            (encodeRowWith (fromJust . packDecimal) . second countEachCell)
          where
            countEachCell :: [BED] -> [(Int, Int)]
            countEachCell = M.toList . foldl' f M.empty
              where
                f m tag = foldl' (\x k -> M.insertWith (+) k (1::Int) x) m $ concat $
                    IM.elems $ IM.containing (M.lookupDefault IM.empty (tag^.chrom) bedTree) p
                  where
                    p | tag^.strand == Just True = tag^.chromStart
                      | tag^.strand == Just False = tag^.chromEnd - 1
                      | otherwise = error "Unkown strand"
        bedTree = bedToTree (++) $ zip regions $ map return [0..]
        header = B.pack $ printf "Sparse matrix: %d x %d" nCell (length regions)

data GeneAccDef = PromoterOnly -- ^ -/+ 1000 around TSS
                | PromoterPlusGeneBody   -- ^ Gene body plus upstream 2000

writePromoters :: SCATACSeqConfig config 
               => GeneAccDef 
               -> ReaderT config IO (File '[] 'Bed)
writePromoters def = do
    dir <- asks ((<> "/Feature/Gene/") . _scatacseq_output_dir) >>= getPath
    let output = dir <> "promoters.bed"
    genes <- asks _scatacseq_annotation >>= liftIO . readGenesValidated . fromJust
    liftIO $ writeBed output $ concatMap getPromoter genes
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
combineExprTable :: SCATACSeqConfig config
                 => FilePath
                 -> ( Maybe (File '[] 'Bed)  -- promoter
                    , [FilePath] )
                 -> ReaderT config IO (Maybe (File '[GeneQuant] 'Tsv))
combineExprTable output (Just promoterFl, fls) = liftIO $ do
    r0 <- decodeFile $ head fls :: IO [(T.Text, U.Vector Int)]
    (clNames, vectors) <- foldM f (unzip r0) $ tail fls
    promoters <- readBed $ promoterFl^.location :: IO [BED]
    let (geneNames, val) = unzip $ map (normalize promoters . U.map fromIntegral) vectors
    DF.writeTable output (T.pack . show) $
        DF.mkDataFrame (head geneNames) clNames $ transpose val
    return $ Just $ location .~ output $ emptyFile
  where
    f (nms, vec) fl = do
        (nms', vec') <- unzip <$> decodeFile fl
        if nms' /= nms
          then error ""
          else return $ force (nms, zipWith (U.zipWith (+)) vec vec')
    normalize promoters vec = unzip $ sort $ M.toList $ M.fromListWith max $
        zipWith g promoters $ U.toList vec
      where
        g bed x = (T.pack $ B.unpack $ fromJust $ bed^.name, x / (fromIntegral (size bed) / 1000) / total)
        total = U.sum vec / 1000000 :: Double
combineExprTable _ _ = return Nothing

-- | Combine expression data into a table and output
mkExprTable :: SCATACSeqConfig config
            => FilePath
            -> ( SCATACSeq S (File '[Gzip] 'Other)  -- matrix
               , File '[] 'Bed  -- promoter
               , SCATACSeq S (File '[] 'Other) )  -- cluster
            -> ReaderT config IO FilePath
mkExprTable prefix (matFl, promoterFl, clFl) = do
    dir <- asks _scatacseq_output_dir >>= getPath . (<> asDir prefix)
    let output = printf "%s/%s_rep%d_gene_count.bin" dir (T.unpack $ matFl^.eid)
            (matFl^.replicates._1)
    liftIO $ do
        let i = B.pack $ T.unpack (matFl^.eid) <> "_" <> show (matFl^.replicates._1)
        promoters <- readBed $ promoterFl^.location
        clusters <- fmap (map (g i)) $ decodeFile $ clFl^.replicates._2.files.location
        mat <- mkSpMatrix readInt $ matFl^.replicates._2.files.location
        let f (nm, xs) = let nm' = B.concat [i, "+", nm] in (nm', xs)
        mkGeneTable promoters (mapRows f mat) clusters >>= encodeFile output
        return output
  where
    g i x = let c = filter ((i `B.isPrefixOf`) . _cell_barcode) $ _cluster_member x
              in x{_cluster_member = c}

{-
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
-}

mkGeneTable :: [BED]   -- ^ Promoter
            -> SpMatrix Int   -- ^ count matrix
            -> [CellCluster]  -- ^ cluster
            -> IO [(T.Text, U.Vector Int)]
mkGeneTable promoters mats clusters = do 
    vectors <- fmap V.fromList $ replicateM (length clusters) $ UM.replicate (length promoters) 0
    let f (bc, counts) = case M.lookup bc bcMap of
            Nothing -> return ()
            Just idx -> do
                let vec = vectors V.! idx
                forM_ counts $ \(i, x) -> UM.modify vec (+x) i
    runResourceT $ runConduit $ streamRows mats .| mapM_C f
    res <- mapM U.unsafeFreeze $ V.toList vectors
    return $ zip (map (T.pack . B.unpack . _cluster_name) clusters) res
  where
    bcMap = M.fromList $ flip concatMap (zip [0..] clusters) $ \(i, cl) ->
        zip (map _cell_barcode $ _cluster_member cl) $ repeat i
{-# INLINE mkGeneTable #-}