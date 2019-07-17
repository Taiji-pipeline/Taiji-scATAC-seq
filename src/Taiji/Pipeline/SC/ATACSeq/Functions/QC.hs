{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE QuasiQuotes #-}
{-# LANGUAGE DataKinds #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.QC
    ( Stat(..)
    , plotStat
    , readStats
    , showStat
    , decodeStat
    , rmAbnoramlFragment
    , rmChrM
    , rmDup
    , baseCoverageStat
    , totalReads
    , tssEnrichment
    ) where

import           Bio.Data.Bam
import           Bio.HTS
import Control.Monad.State.Strict
import Language.Javascript.JMacro
import Bio.Utils.Functions (slideAverage)
import           Bio.Data.Bed
import           Bio.Data.Bed.Types
import Bio.Utils.Misc (readInt, readDouble)
import qualified Data.ByteString.Char8 as B
import qualified Data.Map.Strict as M
import qualified Data.IntMap.Strict as I
import qualified Data.Text as T
import qualified Data.IntervalMap.Strict      as IM
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Unboxed.Mutable as UM

import Taiji.Prelude
import Taiji.Pipeline.SC.ATACSeq.Types
import Taiji.Utils.Plot
import Taiji.Utils.Plot.Vega

data Stat = Stat
    { _cell_barcode :: B.ByteString
    , _dup_rate :: Double
    , _mito_rate :: Double
    , _te :: Double
    , _uniq_reads :: Int }

plotStat :: SCATACSeqConfig config
         => SCATACSeq S (a, b, File '[] 'Tsv )
         -> ReaderT config IO ()
plotStat input = do
    dir <- qcDir
    let output = dir <> "qc_" <> T.unpack (input^.eid) <> "_rep" <>
            show (input^.replicates._1) <> ".html"
    liftIO $ do
        stats <- fmap (filter ((>=100) . _uniq_reads)) $ readStats $
            input^.replicates._2.files._3.location
        let plt = contour $ zip (map (fromIntegral . _uniq_reads) stats) $
                map _te stats
            axes = option [jmacroE| {
                axes: [
                    {
                        domain: false, grid: true, orient: "bottom",
                        scale: "x", title: "number of reads",
                        labelSeparation: 15, labelOverlap: true
                    }, {
                        domain: false, grid: true, orient: "left",
                        scale: "y", title: "TSS enrichment"
                    }
                ]
            } |]
            scales = option [jmacroE| {
                scales: [ {
                    name: "x",
                    type: "log",
                    round: true,
                    nice: true,
                    zero: false,
                    domain: {data: "source", field: "x"},
                    range: "width"
                }, {
                    name: "y",
                    type: "symlog",
                    round: true,
                    nice: true,
                    zero: false,
                    domain: {data: "source", field: "y"},
                    range: "height"
                }, {
                    name: "color",
                    type: "linear",
                    zero: true,
                    domain: {data: "contours", field: "value"},
                    range: "heatmap"
                } ]
            } |]
            hline = option [jmacroE| {
                marks: {
                    type: "rule",
                    encode: {
                        enter: {
                            x2: {signal: "width"},
                            y: {value: 5, scale: "y"},
                            strokeDash: {value: [4,4]}
                        }
                    }
                }
            } |]
            vline = option [jmacroE| {
                marks: {
                    type: "rule",
                    encode: {
                        enter: {
                            y2: {signal: "height"},
                            x: {value: 1000, scale: "x"},
                            strokeDash: {value: [4,4]}
                        }
                    }
                }
            } |]
            n = length $ filter (\x -> _te x >= 5 && _uniq_reads x >= 1000) stats
        savePlots output [plt <> title ("pass QC: " <> show n) <> axes <> vline <> hline <> scales] []

readStats :: FilePath -> IO [Stat]
readStats = fmap (map decodeStat . B.lines) . B.readFile

showStat :: Stat -> B.ByteString
showStat Stat{..} = B.intercalate "\t" $
    [ _cell_barcode
    , B.pack $ show _dup_rate
    , B.pack $ show _mito_rate
    , B.pack $ show _te
    , B.pack $ show _uniq_reads ]
{-# INLINE showStat #-}

decodeStat :: B.ByteString -> Stat
decodeStat x = Stat bc (readDouble dupRate) (readDouble mitoRate)
    (readDouble te) (readInt uniq)
  where
    [bc, dupRate, mitoRate, te, uniq] = B.split '\t' x
{-# INLINE decodeStat #-}

baseCoverageStat :: BAMHeader -> [BAM] -> [(Int, Int)]
baseCoverageStat header bam = sort $ I.toList $ runIdentity $ runConduit $
    mergeBedWith counts beds .| foldlC (I.unionWith (+)) I.empty
  where
    beds = runIdentity $ runConduit $
        yieldMany (sortBy (comparing queryName) bam) .|
        toBed header .| sinkList
    counts xs = I.fromListWith (+) $ flip zip (repeat 1) $
        flip map [lo..hi-1] $ \i -> length $
        filter (\x -> i >= x^.chromStart && i < x^.chromEnd) xs
      where
        lo = minimum $ map (^.chromStart) xs
        hi = maximum $ map (^.chromEnd) xs
{-# INLINE baseCoverageStat #-}

toBed :: Monad m => BAMHeader -> ConduitT BAM BED3 m ()
toBed header = sortedBamToBedPE header .| concatMapC makeFragment
  where
    makeFragment (b1, b2)
        | b1^.chrom /= b2^.chrom || b1^.strand == b2^.strand = Nothing
        | otherwise = 
            let left = if fromJust (b1^.strand)
                    then b1^.chromStart else b2^.chromStart
                right = if not (fromJust $ b2^.strand)
                    then b2^.chromEnd else b1^.chromEnd
            in if left < right
                  then Just (asBed (b1^.chrom) left right :: BED3)
                  else error "Left coordinate is larger than right coordinate."
{-# INLINE toBed #-}

rmAbnoramlFragment:: [BAM] -> State Stat [BAM]
rmAbnoramlFragment = return . filter f  
  where
    f bam = let s = abs $ tLen bam
            in s > 30 && s < 1000
{-# INLINE rmAbnoramlFragment #-}

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
rmChrM :: BAMHeader -> [BAM] -> State Stat ([BAM], [BAM])
rmChrM header input = do
    modify' $ \x -> x{_mito_rate = chrMRate}
    return output
  where
    output = partition notChrM input
    chrMRate = fromIntegral (length $ snd output) / fromIntegral (length input)
    notChrM x = let chr = fromJust $ refName header x
                in chr /= "chrM" && chr /= "M"
{-# INLINE rmChrM #-}

totalReads :: [a] -> State Stat ()
totalReads x = modify' $ \s -> s{_uniq_reads=length x}
{-# INLINE totalReads #-}

tssEnrichment :: BEDTree (Int, Bool)   -- ^ TSS
              -> BAMHeader -> [BAM] -> State Stat ()
tssEnrichment regions header input = modify' $ \x -> x{_te = te}
  where
    te = U.maximum $ normalize vec 
      where
        vec = U.create $ do
            v <- UM.replicate 4000 0
            forM_ input $ \r -> do
                let cutsite = getCutSite r
                forM_ (IM.elems $ intersecting regions cutsite) $ \(x, str) ->
                    UM.unsafeModify v (+1) $ if str
                        then cutsite^.chromStart - (x - 2000)
                        else 3999 - (x + 2000 - cutsite^.chromStart)
            return v
    normalize vec = slideAverage 25 $ U.map (/bk) vec
      where
        bk = (U.sum (U.take 100 vec) + U.sum (U.drop 3900 vec)) / 200 + 0.1
    getCutSite bam = BED3 (fromJust $ refName header bam) i $ i + 1
      where
        i = if isRev bam then endLoc bam - 76 else startLoc bam + 75