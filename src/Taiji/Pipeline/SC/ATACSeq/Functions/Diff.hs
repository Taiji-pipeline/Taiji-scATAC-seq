{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE DataKinds #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Diff
    ( sampleCells
    , mkRefMat
    , diffPeaks
    , diffGenes
    , rpkmPeak
    , rpkmDiffPeak
    ) where

import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U
import qualified Data.ByteString.Char8 as B
import qualified Data.Text as T
import AI.Clustering.Hierarchical
import Bio.Data.Bed
import Data.Binary (decodeFile)
import Bio.Data.Bed.Utils (rpkmBed)
import Data.Conduit.Internal (zipSources)
import System.Random.MWC (create)
import System.Random.MWC.Distributions (uniformShuffle)
import qualified Data.HashSet as S
import qualified Data.HashMap.Strict as M
import qualified Data.Matrix.Unboxed as MU
import Data.Conduit.Zlib (gzip)
import Shelly hiding (FilePath)

import Taiji.Prelude
import Taiji.Pipeline.SC.ATACSeq.Types
import Taiji.Pipeline.SC.ATACSeq.Functions.Utils

{-
sampleCells :: Int   -- ^ number of cells
            -> SCATACSeq S (File tags 'Other)
            -> IO [B.ByteString]
sampleCells n input = do
    mat <- mkSpMatrix id $ input^.replicates._2.files.location
    v <- runResourceT $ runConduit $ streamRows mat .| mapC fst .| sinkVector 
    g <- create
    map f . V.toList . V.take n <$> uniformShuffle v g
  where
    f x = B.pack (T.unpack $ input^.eid) <> "+" <> x
-}

sampleCells :: Int   -- ^ number of cells
            -> [SCATACSeq S (File tags 'Other)]
            -> IO [[B.ByteString]]
sampleCells n input = do
    cls <- decodeFile $ head input^.replicates._2.files.location
    g <- create
    forM cls $ \cl -> do
        let v = V.fromList $ map _cell_barcode $ _cluster_member cl
        V.toList . V.take n <$> uniformShuffle v g

mkRefMat :: SCATACSeqConfig config
         => FilePath
         -> Bool
         -> ( [[B.ByteString]]
            , [SCATACSeq S (File '[Gzip] 'Other)] )
         -> ReaderT config IO FilePath
mkRefMat filename addName (bcs, fls) = do
    dir <- asks _scatacseq_output_dir >>= getPath
    let output = dir <> filename
    liftIO $ do
        mat <- mkSpMatrix id $ head fls^.replicates._2.files.location
        let header = B.pack $ printf "Sparse matrix: %d x %d" (S.size bc) (_num_col mat)
        runResourceT $ runConduit $ mapM_ sourceMat fls .|
            (yield header >> mapC (encodeRowWith id)) .| unlinesAsciiC .|
            gzip .| sinkFile output
        return output
  where
    bc = S.fromList $ concat bcs
    sourceMat input = do
        let filterFn x | addName = (B.pack (T.unpack $ input^.eid) <> "+" <> x) `S.member` bc
                       | otherwise = x `S.member` bc
        mat <- liftIO $ mkSpMatrix id $ input^.replicates._2.files.location
        streamRows mat .| filterC (filterFn . fst)

diffPeaks :: SCATACSeqConfig config
          => ( SCATACSeq S (File tags 'Other)
             , File '[Gzip] 'NarrowPeak
             , FilePath )
          -> ReaderT config IO (SCATACSeq S (File '[Gzip] 'NarrowPeak))
diffPeaks (input, peakFl, ref) = do
    dir <- asks ((<> "/Diff/Peak/") . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d.np.gz" dir
            (T.unpack $ input^.eid) (input^.replicates._1)
    input & replicates.traversed.files %%~ liftIO . ( \fl -> do
        stats <- fmap (M.fromList . map (\(a,b,c,d) -> (a, (b,c,d))) . filter (\x -> x^._4 < 0.01)) $
            diffAnalysis (fl^.location) ref
        let f (i, peak) = case M.lookup i stats of
                Nothing -> Nothing
                Just (fold, pval, fdr) -> 
                    let logP | pval == 0 = 200
                             | otherwise = negate $ logBase 10 pval 
                        logFDR | fdr == 0 = 200
                               | otherwise = negate $ logBase 10 fdr
                    in Just $ npSignal .~ fold $
                        npPvalue .~ Just logP $
                        npQvalue .~ Just logFDR $ peak
        runResourceT $ runConduit $
            zipSources (iterateC succ 0) (streamBedGzip $ peakFl^.location) .|
            concatMapC f .| sinkFileBedGzip output
        return $ location .~ output $ emptyFile )


-- | Compute RPKM for each peak
rpkmPeak :: SCATACSeqConfig config
         => ( (B.ByteString, File '[Gzip] 'Bed)
            , Maybe (File '[Gzip] 'NarrowPeak) )
         -> ReaderT config IO (B.ByteString, FilePath)
rpkmPeak ((nm, tags), Just peakFl) = do
    dir <- asks _scatacseq_output_dir >>= getPath . (<> "/Feature/Peak/")
    let output = dir <> "peak_rpkm_" <> B.unpack nm <> ".txt"
    liftIO $ do
        peaks <- runResourceT $ runConduit $
            streamBedGzip (peakFl^.location) .| sinkList :: IO [BED3]
        vec <- runResourceT $ runConduit $ streamBedGzip (tags^.location) .| rpkmBed peaks 
        B.writeFile output $ B.unlines $ map toShortest $ U.toList vec
        return (nm, output)
rpkmPeak _ = undefined

rpkmDiffPeak :: SCATACSeqConfig config
             => ( [SCATACSeq S (File '[Gzip] 'NarrowPeak)]   -- ^ Diff peaks
                , Maybe (File '[Gzip] 'NarrowPeak)  -- ^ All peaks
                , [(B.ByteString, FilePath)] )  -- ^ RPKM
             -> ReaderT config IO FilePath
rpkmDiffPeak (diffs, Just peak, rpkms) = do
    dir <- asks ((<> "/Diff/Peak/") . _scatacseq_output_dir) >>= getPath
    let output = dir <> "diff_peak_RPKM.tsv"
    liftIO $ do
        mat <- readRPKMs $ map snd rpkms
        peakList <- runResourceT $ runConduit $
            streamBedGzip (peak^.location) .| sinkList :: IO [BED3]
        res <- forM diffs $ \input -> do
            peaks <- fmap (filter f) $ runResourceT $ runConduit $
                streamBedGzip (input^.replicates._2.files.location) .| sinkList
            let idx = getPeakIndex peaks peakList
                submat = map (B.intercalate "\t" . map toShortest . U.toList) $
                    MU.toRows $ getSubMatrix idx mat
            return $ "#" <> B.pack (T.unpack $ input^.eid) : submat
        B.writeFile output $ B.unlines $ header : concat res
        return output
  where
    f x = x^.npSignal > 1.5 || x^.npSignal < 1/1.5
    header = B.intercalate "\t" $ fst $ unzip rpkms

getSubMatrix :: [Int] -> MU.Matrix Double -> MU.Matrix Double
getSubMatrix idx mat = MU.fromRows $ flatten $ hclust Ward dat euclidean
  where
    dat = V.fromList $ map (mat `MU.takeRow`) idx

-- | Get the indices of query peaks in the reference peak list.
getPeakIndex :: (BEDLike b1, BEDLike b2)
             => [b1]   -- ^ Differential peaks
             -> [b2]   -- ^ Reference peaks
             -> [Int]  -- ^ Indices
getPeakIndex query ref = flip map query $ \q -> M.lookupDefault undefined
    (q^.chrom, q^.chromStart, q^.chromEnd) peakIdxMap
  where
    peakIdxMap = M.fromList $ zipWith
        (\x i -> ((x^.chrom, x^.chromStart, x^.chromEnd), i)) ref [0..]

readRPKMs :: [FilePath] -> IO (MU.Matrix Double)
readRPKMs fls = do
    vecs <- mapM readData fls
    return $ MU.map (logBase 2 . (+1)) $ MU.fromColumns vecs
  where
    readData fl = U.fromList . map readDouble . B.lines <$> B.readFile fl

{-
mkDiffPeakFig :: SCATACSeqConfig config 
              => FilePath
              -> ReaderT config IO ()
mkDiffPeakFig fl = do
    dir <- figDir
    (header:rest) <- B.lines <$> B.readFile fl
-}


diffGenes :: SCATACSeqConfig config
          => ( SCATACSeq S (File tags 'Other)
             , FilePath   -- ^ Gene names
             , FilePath ) -- ^ Ref matrix
          -> ReaderT config IO (SCATACSeq S (File '[] 'Tsv))
diffGenes (input, nameFl, ref) = do
    dir <- asks ((<> "/Diff/Gene/") . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d.tsv" dir
            (T.unpack $ input^.eid) (input^.replicates._1)
    input & replicates.traversed.files %%~ liftIO . ( \fl -> do
        stats <- fmap (M.fromList . map (\(a,b,c,d) -> (a, (b,c,d))) . filter (\x -> x^._4 < 0.01)) $
            diffAnalysis (fl^.location) ref
        let f (i, nm) = case M.lookup i stats of
                Nothing -> Nothing
                Just (fold, pval, fdr) ->
                    let logP | pval == 0 = 200
                             | otherwise = negate $ logBase 10 pval 
                        logFDR | fdr == 0 = 200
                               | otherwise = negate $ logBase 10 fdr
                    in Just $ B.intercalate "\t" [nm, toShortest fold,
                            toShortest logP, toShortest logFDR]
        runResourceT $ runConduit $ zipSources (iterateC succ 0)
            (sourceFile nameFl .| linesUnboundedAsciiC .| mapC (head . B.words)) .|
            concatMapC f .| unlinesAsciiC .| sinkFile output
        return $ location .~ output $ emptyFile )


diffAnalysis :: FilePath  -- ^ foreground
             -> FilePath  -- ^ background
             -> IO [(Int, Double, Double, Double)]
diffAnalysis fg bk = withTemp Nothing $ \tmp -> do
    shelly $ run_ "sc_utils" [ "diff", T.pack tmp, "--fg", T.pack fg
        , "--bg", T.pack bk ]
    map (f . B.words) . B.lines <$> B.readFile tmp
  where
    f [a,b,c,d] = (readInt a, readDouble b, readDouble c, readDouble d)
    f _ = error "wrong format!"

