{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE RecordWildCards #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Align
    ( tagAlign
    , filterBamSort
    , rmDuplicates
    ) where

import Bio.Data.Experiment
import Bio.Data.Bed
import Bio.Data.Bam
import           Bio.HTS
import qualified Data.ByteString.Char8 as B
import Conduit
import Data.Conduit.Internal (zipSinks)
import Control.Monad.State.Strict
import Control.Lens
import Data.Conduit.List (groupBy)
import Bio.Pipeline hiding (frip)
import Data.Function (on)
import Data.Either (fromRight)
import qualified Data.Text as T
import qualified Data.Map.Strict as M
import Data.Maybe
import Text.Printf (printf)
import Scientific.Workflow
import Control.Monad.Reader (asks, liftIO)

import Taiji.Pipeline.SC.ATACSeq.Types
import Taiji.Pipeline.SC.ATACSeq.Functions.Utils

tagAlign :: SCATACSeqConfig config
         => SCATACSeq S (SomeTags 'Fastq, SomeTags 'Fastq)
         -> WorkflowConfig config ( SCATACSeq S (File '[PairedEnd] 'Bam) )
tagAlign input = do
    dir <- asks ((<> "/Bam") . _scatacseq_output_dir) >>= getPath
    idx <- asks (fromJust . _scatacseq_bwa_index)
    let output = printf "%s/%s_rep%d.bam" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . ( \(f1,f2) ->
        let f1' = fromSomeTags f1 :: File '[] 'Fastq
            f2' = fromSomeTags f2 :: File '[] 'Fastq
        in fromRight undefined <$> bwaAlign output idx (Right (f1', f2'))
            (defaultBWAOpts & bwaCores .~ 8)
        )

filterBamSort :: SCATACSeqConfig config
              => SCATACSeq S (File '[PairedEnd] 'Bam)
              -> WorkflowConfig config
                  (SCATACSeq S (File '[NameSorted, PairedEnd] 'Bam))
filterBamSort input = do
    dir <- asks ((<> "/Bam") . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_filt.bam" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . filterBam "./" output

rmDuplicates :: SCATACSeqConfig config
             => SCATACSeq S (File '[NameSorted, PairedEnd] 'Bam)
             -> WorkflowConfig config
                (SCATACSeq S (File '[NameSorted, Gzip] 'Bed, File '[] 'Tsv))
rmDuplicates input = do
    dir <- asks _scatacseq_output_dir >>= getPath . (<> (asDir "/Bed"))
    anno <- fromJust <$> asks _scatacseq_annotation 
    let output1 = printf "%s/%s_rep%d_filt.bed.gz" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
        output2 = printf "%s/%s_rep%d_stat.txt" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
        f fl = do
            header <- getBamHeader fl
            promoters <- (fmap.fmap) (const ()) <$> readPromoters anno
            _ <- runResourceT $ runConduit $ streamBam fl .|
                groupBy ((==) `on` (extractBarcode . queryName)) .|
                mapC (filterReads header promoters) .|
                zipSinks 
                    (concatMapC filterCell .| concatC .| sinkFileBedGzip output1)
                    (mapC (showStat . snd) .| unlinesAsciiC .| sinkFile output2)
            return ( location .~ output1 $ emptyFile
                   , location .~ output2 $ emptyFile )
    input & replicates.traverse.files %%~ liftIO . f . (^.location)

filterCell :: (a, Stat) -> Maybe a 
filterCell (x, Stat{..})
    | _uniq_reads >= 1000 && _frip <= 0.8 && _frip >= 0.2 = Just x
    | otherwise = Nothing
{-# INLINE filterCell #-}

data Stat = Stat
    { _cell_barcode :: B.ByteString
    , _dup_rate :: Double
    , _mito_rate :: Double
    , _frip :: Double
    , _uniq_reads :: Int }

showStat :: Stat -> B.ByteString
showStat (Stat a b c d e) = B.intercalate "\t" $
    [ a, B.pack $ show b, B.pack $ show c, B.pack $ show d, B.pack $ show e]
{-# INLINE showStat #-}

-- | Remove duplicated and chrM reads. 
filterReads :: BAMHeader
            -> BEDTree a   -- ^ Promoter
            -> [BAM]
            -> ([BED], Stat)
filterReads hdr promoters bam = runState f $ Stat bc 0 0 0 0
  where
    f = mapMaybe (fmap changeName . bamToBed hdr) <$> rmDup bam >>=
        rmChrM >>= frip promoters >>= final
    final x = do
        modify' $ \s -> s{_uniq_reads=length x}
        return x
    bc = extractBarcode $ queryName $ head bam
    changeName x = name .~ Just bc $ x
{-# INLINE filterReads #-}

-- | Remove duplicates for reads originated from a single cell.
rmDup :: [BAM] -> State Stat [BAM]
rmDup input = do
    modify' $ \x -> x{_dup_rate = dupRate}
    return output
  where
    output = map snd $ M.elems $ M.fromListWith collapse $
        map (\b -> (getKey b, (getSc b, b))) input
    collapse x y | fst x > fst y = x
                 | otherwise = y
    getKey b | not (hasMultiSegments flg) || isNextUnmapped flg = singleKey
             | otherwise = pairKey
      where
        (singleKey, pairKey) = makeKey (const Nothing) b
        flg = flag b
    getSc b = case queryAuxData ('m', 's') b of
        Just (AuxInt x) -> x + fromJust (sumQual 15 b)
        _ -> fromJust $ sumQual 15 b
    dupRate = 1 - fromIntegral (length output) / fromIntegral (length input)
{-# INLINE rmDup #-}

-- | Remove chrM reads.
rmChrM :: [BED] -> State Stat [BED]
rmChrM input = do
    modify' $ \x -> x{_mito_rate = chrMRate}
    return output
  where
    output = filter notChrM input
    chrMRate = 1 - fromIntegral (length output) / fromIntegral (length input)
    notChrM x = x^.chrom /= "chrM" && x^.chrom /= "M"
{-# INLINE rmChrM #-}

-- | Fraction of Reads in Promoter
frip :: BEDTree a -> [BED] -> State Stat [BED]
frip promoter input = do
    modify' $ \x -> x{_frip=rate}
    return input
  where
    rate = fromIntegral (length $ filter (isIntersected promoter) input) /
        fromIntegral (length input)
{-# INLINE frip #-}
