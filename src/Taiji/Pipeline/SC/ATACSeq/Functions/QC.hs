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
    , readTSS
    , detectDoublet
    , removeDoublet
    , plotClusterQC
    ) where

import           Bio.Data.Bam
import           Bio.HTS
import Control.Monad.State.Strict
import Language.Javascript.JMacro
import Bio.Utils.Functions (slideAverage)
import Data.Binary (decodeFile)
import           Bio.Data.Bed
import Bio.RealWorld.GENCODE
import           Bio.Data.Bed.Types
import Bio.Utils.Misc (readInt, readDouble)
import Control.Arrow (second)
import qualified Data.ByteString.Char8 as B
import qualified Data.Map.Strict as M
import qualified Data.IntMap.Strict as I
import qualified Data.HashSet as S
import Data.Conduit.List (groupBy)
import Data.Conduit.Internal (zipSinks)
import qualified Data.Text as T
import Data.List.Ordered (nubSort)
import qualified Data.IntervalMap.Strict      as IM
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Unboxed.Mutable as UM
import Shelly (shelly, run_)

import Taiji.Prelude hiding (groupBy)
import Taiji.Pipeline.SC.ATACSeq.Types
import Taiji.Utils.Plot
import Taiji.Utils.Plot.Vega
import qualified Taiji.Utils.Plot.ECharts as E
import Taiji.Pipeline.SC.ATACSeq.Functions.Utils (mkSpMatrix, streamRows)
import qualified Taiji.Utils.DataFrame as DF

data Stat = Stat
    { _barcode :: B.ByteString
    , _dup_rate :: Double
    , _mito_rate :: Double
    , _te :: Double
    , _uniq_reads :: Int
    , _doublet_score :: Double }

plotStat :: SCATACSeqConfig config
         => [SCATACSeq S (a, b, File '[] 'Tsv )]
         -> ReaderT config IO ()
plotStat inputs = do
    dir <- qcDir
    let output = dir <> "/qc.html"
    liftIO $ do
        (cellQCs, stats) <- fmap unzip $ forM inputs $ \input -> do
            stats <- readStats $ input^.replicates._2.files._3.location
            let stats' = filter (\x -> _te x >= 7 && _uniq_reads x >= 1000) stats
                cellQC = plotCells stats <> title
                    (printf "%s: %d cells passed QC" (T.unpack $ input^.eid) (length stats'))
            return (cellQC, (input^.eid, stats'))
        savePlots output ([plotNumReads stats, plotDupRate stats, plotMitoRate stats] <> cellQCs) []

readStats :: FilePath -> IO [Stat]
readStats = fmap (map decodeStat . B.lines) . B.readFile

showStat :: Stat -> B.ByteString
showStat Stat{..} = B.intercalate "\t" $
    [ _barcode
    , B.pack $ show _dup_rate
    , B.pack $ show _mito_rate
    , B.pack $ show _te
    , B.pack $ show _uniq_reads
    , B.pack $ show _doublet_score ]
{-# INLINE showStat #-}

decodeStat :: B.ByteString -> Stat
decodeStat x = Stat bc (readDouble dupRate) (readDouble mitoRate)
    (readDouble te) (readInt uniq) (readDouble ds)
  where
    [bc, dupRate, mitoRate, te, uniq, ds] = B.split '\t' x
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

totalReads :: [BAM] -> State Stat ()
totalReads x = modify' $ \s -> s
    {_uniq_reads = length $ filter (isFirstSegment . flag) x}
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
    normalize vec = slideAverage 5 $ U.map (/bk) vec
      where
        bk = (U.sum (U.take 100 vec) + U.sum (U.drop 3900 vec)) / 200 + 0.1
    getCutSite bam = BED3 (fromJust $ refName header bam) i $ i + 1
      where
        i = if isRev bam then endLoc bam - 76 else startLoc bam + 75

readTSS :: FilePath -> IO (BEDTree (Int, Bool))
readTSS = fmap (bedToTree const . concatMap fn) . readGenes
  where
    fn Gene{..} = map g $ nubSort tss
      where
        g x = (BED3 geneChrom (x - 2000) (x + 2000), (x, geneStrand))
        tss | geneStrand = geneLeft : map fst geneTranscripts
            | otherwise = geneRight : map snd geneTranscripts

removeDoublet :: SCATACSeqConfig config
              => SCATACSeq S ( File '[NameSorted, Gzip] 'Bed
                             , (File '[] 'Tsv, Double) )
              -> ReaderT config IO (SCATACSeq S (File '[NameSorted, Gzip] 'Bed, Int))
removeDoublet input = do
    dir <- asks _scatacseq_output_dir >>= getPath . (<> (asDir "/Bed"))
    let output = printf "%s/%s_rep%d_srt_filt_no_dblet.bed.gz" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . (\(bedFl, (statFl, th)) -> do
        stats <- readStats $ statFl^.location
        let cells = S.fromList $ map _barcode $ filter ((<=th) . _doublet_score) stats
            f :: [BED] -> Bool
            f = (`S.member` cells) . fromJust . (^.name) . head
            sink = zipSinks lengthC $ concatC .| sinkFileBedGzip output
        (nCells, _) <- runResourceT $ runConduit $ streamBedGzip (bedFl^.location) .|
            groupBy ((==) `on` (^.name)) .| filterC f .| sink
        return (location .~ output $ emptyFile, nCells) )

detectDoublet :: SCATACSeqConfig config
              => SCATACSeq S (File tags 'Other, (a, b, File '[] 'Tsv ))
              -> ReaderT config IO (SCATACSeq S (File '[] 'Tsv, Double))
detectDoublet input = do
    dir <- qcDir
    let output = dir <> "doublet_" <> T.unpack (input^.eid) <> "_rep" <>
            show (input^.replicates._1) <> ".txt"
        outputPlot = dir <> "doublet_" <> T.unpack (input^.eid) <> "_rep" <>
            show (input^.replicates._1) <> ".html"
    input & replicates.traverse.files %%~ liftIO . (\(fl, (_,_,stat)) -> do
        shelly $ run_ "sc_utils" ["doublet", T.pack $ fl^.location, T.pack output]
        [thres, sc, sim_sc] <- B.lines <$> B.readFile output
        let th = read $ B.unpack thres :: Double
            ds = map readDouble $ B.words sc
            ds_sim = map readDouble $ B.words sim_sc
            rate = fromIntegral (length $ filter (>=th) ds) /
                fromIntegral (length ds) * 100 :: Double
        savePlots outputPlot [mkHist ds th <> title (printf "doublet percentage: %.1f%%" rate)
            , mkHist ds_sim th] []

        statMap <- fmap (M.fromList . map (\x -> (_barcode x, x))) $ readStats $ stat^.location
        mat <- mkSpMatrix id $ fl^.location
        bcs <- runResourceT $ runConduit $ streamRows mat .| mapC fst .| sinkList
        B.writeFile (stat^.location) $ B.unlines $ flip map (zip bcs ds) $ \(bc, val) -> 
            let s = M.findWithDefault undefined bc statMap
            in showStat s{_doublet_score = val}
        return (stat, th)
        )
  where
    mkHist xs ref = plt <> rule
      where
        plt = hist xs 100
        rule = option [jmacroE| {
            layer: {
                mark: "rule",
                data: {values: [{ref: `ref`}]},
                encoding: {
                    x: { field: "ref"},
                    color: {value: "red"},
                    size: {value: 1}
                }
            }
       } |]


--------------------------------------------------------------------------------
-- Plotting
--------------------------------------------------------------------------------

--plotNumCells :: [(T.Text, [Stat])] -> Vega
--plotNumCells stats = $ flip map stats $ \(nm, stat) -> (nm, length stat)

plotDupRate :: [(T.Text, [Stat])] -> Vega
plotDupRate stats = violin $ flip map stats $ \(nm, stat) -> 
    (nm, map _dup_rate stat)

plotMitoRate :: [(T.Text, [Stat])] -> Vega
plotMitoRate stats = violin $ flip map stats $ \(nm, stat) -> 
    (nm, map _mito_rate stat)

plotNumReads :: [(T.Text, [Stat])] -> Vega
plotNumReads stats = violin $ flip map stats $ \(nm, stat) -> 
    (nm, map (logBase 10 . fromIntegral . _uniq_reads) stat)

plotCells :: [Stat] -> Vega
plotCells input = plt <> axes <> vline <> hline <> scales
  where
    plt = contour $ zip (map (fromIntegral . _uniq_reads) stats) $ map _te stats
    stats = filter ((>=100) . _uniq_reads) input
    axes = option [jmacroE| {
        axes: [
            {
                domain: false, grid: true, orient: "bottom",
                scale: "x", title: "number of fragments",
                labelSeparation: 15, labelOverlap: true
            }, {
                domain: false, grid: true, orient: "left",
                scale: "y", title: "TSS enrichment",
                labelSeparation: 5, labelOverlap: true
            }
        ]
    } |]
    scales = option [jmacroE| {
        scales: [ {
            name: "x", type: "log", round: true, nice: true,
            zero: false, domain: {data: "source", field: "x"}, range: "width"
        }, {
            name: "y", type: "linear", round: true, nice: true, zero: false,
            domain: {data: "source", field: "y"}, range: "height"
        }, {
            name: "color", type: "linear", zero: true,
            domain: {data: "contours", field: "value"}, range: "heatmap"
        } ]
    } |]
    hline = option [jmacroE| {
        marks: {
            type: "rule",
            encode: {
                enter: {
                    x2: {signal: "width"},
                    y: {value: 7, scale: "y"},
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


plotClusterQC :: SCATACSeqConfig config
              => SCATACSeq S
                ( (File '[] 'Tsv, a)
                , File '[] 'Other
                , (FilePath, File '[Gzip] 'Other)
                , ([String], File '[] 'Tsv) )
              -> ReaderT config IO ()
plotClusterQC input = do
    dir <- qcDir
    let output = dir <> T.unpack (input^.eid) <> "_cluster_qc.html"
    liftIO $ do
        stats <- readStats $ statFl^.location
        cls <- decodeFile $ clFl^.location
        geneExpr <- M.fromList <$> readGeneExpr (map B.pack genes) idxFl matFl
        df <- DF.readTable $ markerFl^.location
        clusterQC output cls stats (genes, geneExpr) df
  where
    ((statFl,_), clFl, (idxFl, matFl), (genes, markerFl)) = input^.replicates._2.files

clusterQC :: FilePath
          -> [CellCluster]
          -> [Stat]
          -> ([String], M.Map B.ByteString [Int])
          -> DF.DataFrame Double
          -> IO ()
clusterQC output cs stats (names, expr) df = savePlots output []
    [ E.scatter3D' dat3D viz <> E.toolbox
    , hp <> E.toolbox <> E.option
        [jmacroE| { visualMap: { inRange: {color: `viridis`} } } |]
    ]
  where
    viz = statViz <> geneViz
    statViz = [ ("read depth", map (logBase 10 . fromIntegral . _uniq_reads) statOrdered)
              , ("TSS enrichment", map _te statOrdered)
              , ("duplication rate", map _dup_rate statOrdered)
              , ("mito rate", map _mito_rate statOrdered)
              , ("doublet score", map _doublet_score statOrdered) ]
    geneViz = zipWith (\x y -> (x, map fromIntegral y)) names $ transpose $
        concatMap (map (\x -> M.findWithDefault undefined (_cell_barcode x) expr) .
            _cluster_member) cs
    dat3D = flip map cs $ \(CellCluster nm cells) ->
        (B.unpack nm, map _cell_3d cells)
    statOrdered = concatMap (map
        (\x -> M.findWithDefault undefined (_cell_barcode x) stats') .
        _cluster_member) cs
    stats' = M.fromList $ map (\x -> (_barcode x, x)) stats

    hp = E.heatmap $ DF.reorderRows (DF.orderByCluster id) $
        DF.reorderColumns (DF.orderByCluster id) df
    viridis :: [String]
    viridis = ["#440154", "#482173", "#433E85", "#38598C", "#2D708E",
        "#25858E", "#1E9B8A", "#2BB07F", "#51C56A", "#85D54A", "#C2DF23", "#FDE725"]
 

readGeneExpr :: [B.ByteString] 
             -> FilePath   -- ^ Index
             -> File '[Gzip] 'Other  -- ^ gene matrix
             -> IO [(B.ByteString, [Int])]
readGeneExpr genes idxFl fl = do
    allGenes <- map (head . B.split '\t') . B.lines <$> B.readFile idxFl
    let idxMap = M.fromList $ zip allGenes [0..]
        idx = map (\x -> M.findWithDefault (error $ show x) x idxMap) genes
        f xs = let m = M.fromList xs
               in map (\k -> M.findWithDefault 0 k m) idx
    mat <- mkSpMatrix readInt $ fl^.location
    runResourceT $ runConduit $ streamRows mat .| mapC (second f) .| sinkList