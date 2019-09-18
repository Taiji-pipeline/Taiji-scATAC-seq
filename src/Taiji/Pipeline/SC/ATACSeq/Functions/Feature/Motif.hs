{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE PartialTypeSignatures #-}

module Taiji.Pipeline.SC.ATACSeq.Functions.Feature.Motif
    ( getOpenRegion
    , findMotifsPre
    , findMotifs
    ) where

import Bio.Seq.IO (withGenome, getChrSizes)
import Bio.Data.Bed.Utils
import           Bio.Pipeline
import           Bio.Pipeline.Instances ()
import Bio.Data.Bed
import Data.Binary
import           Bio.Motif                     hiding (score)
import qualified Data.HashMap.Strict as M
import qualified Data.ByteString.Char8 as B
import qualified Bio.Utils.BitVector as BV
import           System.IO (hClose)
import           System.IO.Temp (withTempFile)
import           Data.Default (def)

import Taiji.Prelude
import Taiji.Pipeline.SC.ATACSeq.Types

findMotifsPre :: SCATACSeqConfig config
              => Double   -- ^ Motif finding cutoff
              -> Maybe (File '[Gzip] 'NarrowPeak)   -- ^ Regions to look for
              -> ReaderT config IO
                  [(B.ByteString, File '[Gzip] 'NarrowPeak, File '[] 'Other)]
findMotifsPre p (Just region) = do
    motifFile <- fromMaybe (error "Motif file is not specified!") <$>
        asks _scatacseq_motif_file
    genome <- getGenomeIndex
    chrs <- liftIO $ withGenome genome $ return . map fst . getChrSizes
    dir <- asks _scatacseq_output_dir >>= getPath . (<> (asDir "/temp"))
    let output = dir ++ "/motif.bin"
    liftIO $ map (mkCutoffMotif def p) <$> readMEME motifFile >>= encodeFile output
    return $ zip3 chrs (repeat region) $ repeat $ location .~ output $ emptyFile
findMotifsPre _ _ = return []

-- | Identify motif binding sites.
findMotifs :: SCATACSeqConfig config
           => (B.ByteString, File '[Gzip] 'NarrowPeak, File '[] 'Other)
           -> ReaderT config IO (B.ByteString, Maybe (File '[] 'BigBed))
findMotifs (chr, openChromatin, motifFl) = do
    seqIndex <- getGenomeIndex
    chrSize <- liftIO $ withGenome seqIndex $
        return . filter ((==chr). fst) . getChrSizes
    dir <- asks _scatacseq_output_dir >>= getPath . (<> (asDir "/Feature/TF/"))
    let output = dir ++ B.unpack chr ++ ".bb"
    liftIO $ withGenome seqIndex $ \g -> withTempFile "./" "tmp_tbfs.bed" $ \tmpFl h -> do
        hClose h
        motifs <- decodeFile (motifFl^.location) :: IO [CutoffMotif]
        r <- runResourceT $ runConduit $
            (streamBedGzip (openChromatin^.location) :: _ _ BED3 _ _) .|
            filterC (\x -> x^.chrom == chr) .| scanMotif g motifs .| sink tmpFl
        fl <- if r
            then return Nothing
            else fmap Just $ bedToBigBed output chrSize $
                location .~ tmpFl $ emptyFile
        return (chr, fl)
  where
    sink fl = do
        r <- nullC
        sinkFileBed fl
        return r
{-# INLINE findMotifs #-}

-- | Identify all accessiable regions.
getOpenRegion :: SCATACSeqConfig config
              => [SCATACSeq S (File '[NameSorted, Gzip] 'Bed)]
              -> ReaderT config IO (Maybe (File '[Gzip] 'Bed))
getOpenRegion inputs 
    | null inputs = return Nothing
    | otherwise = do
        dir <- asks _scatacseq_output_dir >>= getPath
        genome <- asks (fromJust . _scatacseq_genome_index)
        chrSize <- liftIO $ withGenome genome $ return . getChrSizes
        let output = dir <> "/accessible_region.bed.gz"
            source = forM_ inputs $ \input ->
                streamBedGzip $ input^.replicates._2.files.location
        runResourceT $ runConduit $ source .|
            (baseMap chrSize >>= baseMapToRegion) .| sinkFileBedGzip output
        return $ Just $ emptyFile & location .~ output
{-# INLINE getOpenRegion #-}

baseMapToRegion :: Monad m => BaseMap -> ConduitT i BED3 m ()
baseMapToRegion (BaseMap bm) = forM_ (M.toList bm) $ \(chr, bv) ->
    mkBedWindow 50 chr bv
{-# INLINE baseMapToRegion #-}

-- | Create bed from basemap.
mkBedWindow :: Monad m 
            => Int   -- ^ half window size
            -> B.ByteString
            -> BV.BitVector -> ConduitT i BED3 m ()
mkBedWindow window chr bv = go (0,0) 0
  where
    go (lo, hi) i
        | i >= n = yield $ asBed chr lo hi
        | not (bv BV.! i) = go (lo, hi) $! i + 1
        | cur_lo <= hi = go (lo, cur_hi) $! i + 1
        | otherwise = yield (asBed chr lo hi) >> go (cur_lo, cur_hi) (i + 1)
      where
        (cur_lo, cur_hi) = getRegion i
    n = BV.size bv
    getRegion x = (max 0 $ x - window, min n $ x + 1 + window)
{-# INLINE mkBedWindow #-}
