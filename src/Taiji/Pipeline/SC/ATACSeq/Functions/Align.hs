{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Align
    ( tagAlign
    , filterBamSort
    , rmDuplicates
    , mkBedGzip
    ) where

import Bio.Data.Experiment
import Bio.Data.Bed
import Bio.Data.Bam
import           Bio.HTS
import qualified Data.ByteString.Char8 as B
import Conduit
import Data.Conduit.Internal (zipSinks)
import Control.Lens
import Data.Conduit.List (groupBy)
import Bio.Pipeline
import Data.Function (on)
import Data.Either (fromRight)
import qualified Data.Text as T
import qualified Data.Map.Strict as M
import Data.Maybe
import Text.Printf (printf)
import Scientific.Workflow
import Control.Monad.Reader (asks, liftIO)

import Taiji.Pipeline.SC.ATACSeq.Types
import Taiji.Pipeline.SC.ATACSeq.Functions.Utils (extractBarcode)

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

mkBedGzip :: SCATACSeqConfig config
          => (SCATACSeq S (File '[NameSorted] 'Bam, File '[] 'Tsv, Int))
          -> WorkflowConfig config
              (SCATACSeq S (File '[NameSorted, Gzip] 'Bed))
mkBedGzip input = do
    dir <- asks ((<> "/Bed") . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d.bed.gz" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ (\(fl,_,_) -> liftIO $ do
        header <- getBamHeader $ fl^.location
        runResourceT $ runConduit $ streamBam (fl^.location) .|
            bamToBedC header .| mapC changeName .|
            sinkFileBedGzip output
        return $ location .~ output $ emptyFile )
  where
    changeName x = name %~ fmap extractBarcode $ x

rmDuplicates :: SCATACSeqConfig config
             => SCATACSeq S (File '[NameSorted, PairedEnd] 'Bam)
             -> WorkflowConfig config
                (SCATACSeq S (File '[NameSorted] 'Bam, File '[] 'Tsv, Int))
rmDuplicates input = do
    dir <- asks _scatacseq_output_dir >>= getPath . (<> (asDir "/Bam"))
    let output1 = printf "%s/%s_rep%d_filt_dedup.bam" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
        output2 = printf "%s/%s_rep%d_stat.txt" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
        f fl = do
            header <- getBamHeader fl
            ((n, _), _) <- runResourceT $ runConduit $ streamBam fl .|
                groupBy ((==) `on` (extractBarcode . queryName)) .|
                mapC (deDupAndRmChrM header) .| zipSinks 
                    (concatMapC filterFn .| zipSinks lengthC (concatC .| sinkBam output1 header))
                    (mapC (toString . snd) .| unlinesAsciiC .| sinkFile output2)
            return ( location .~ output1 $ emptyFile
                   , location .~ output2 $ emptyFile
                   , n )
    input & replicates.traverse.files %%~ liftIO . f . (^.location)
  where
    toString (a,b,c,d) = B.intercalate "\t" $ a : map (B.pack . show) [b,c,d]
    filterFn (x, (_,_,_,n)) | n >= 1000 = Just x
                            | otherwise = Nothing

-- | Remove duplicated and chrM reads. 
deDupAndRmChrM :: BAMHeader
               -> [BAM]
               -> ([BAM], (B.ByteString , Double, Double, Double))
deDupAndRmChrM hdr bam =
    (noChrM, (extractBarcode $ queryName $ head bam, dupRate, chrMRate, nRemain))
  where
    dedup = rmDup bam
    nDedup = fromIntegral $ length dedup
    dupRate = 1 - nDedup / n
    noChrM = rmChrM hdr dedup
    nRemain = fromIntegral $ length noChrM
    chrMRate = 1 - nRemain / nDedup
    n = fromIntegral $ length bam

-- | Remove duplicates for reads originated from a single cell.
rmDup :: [BAM] -> [BAM]
rmDup bam = map snd $ M.elems $ M.fromListWith collapse $
    map (\b -> (getKey b, (getSc b, b))) bam
  where
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
{-# INLINE rmDup #-}

-- | Remove chrM reads.
rmChrM :: BAMHeader -> [BAM] -> [BAM]
rmChrM hdr = filter $ \x ->
    let chr = fromJust $ refName hdr x in chr /= "chrM" && chr /= "M"
{-# INLINE rmChrM #-}
