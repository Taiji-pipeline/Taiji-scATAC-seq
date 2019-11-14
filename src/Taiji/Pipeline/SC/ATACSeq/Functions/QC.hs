{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE QuasiQuotes #-}
{-# LANGUAGE DataKinds #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.QC
    ( Stat(..)
    , passedQC
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
    , plotNumReads
    , plotDupRate
    , plotMitoRate
    , plotTE
    ) where

import           Bio.Data.Bam
import           Bio.HTS
import Control.Monad.State.Strict
import Language.Javascript.JMacro
import Bio.Utils.Functions (slideAverage)
import           Bio.Data.Bed
import Bio.RealWorld.GENCODE
import           Bio.Data.Bed.Types
import Bio.Utils.Misc (readInt, readDouble)
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
import Taiji.Utils (mkSpMatrix, streamRows)

data Stat = Stat
    { _barcode :: B.ByteString
    , _dup_rate :: Double
    , _mito_rate :: Double
    , _te :: Double
    , _uniq_reads :: Int
    , _doublet_score :: Double }

-- | Whether the given cell passes the QC.
passedQC :: Double -> Stat -> Bool
passedQC cutoff x = _te x >= cutoff && _uniq_reads x >= 1000 && _doublet_score x <= 0.5

plotStat :: SCATACSeqConfig config
         => [SCATACSeq S (File '[] 'Tsv)]
         -> ReaderT config IO FilePath
plotStat [] = return ""
plotStat inputs = do
    dir <- qcDir
    teCutoff <- asks _scatacseq_te_cutoff
    let output = dir <> "/qc.html"
        outputStat = dir <> "/qc_stats.tsv"
    liftIO $ do
        (cellQCs, stats) <- fmap unzip $ forM inputs $ \input -> do
            stats <- readStats $ input^.replicates._2.files.location
            let stats' = filter (passedQC teCutoff) stats
                cellQC = plotCells teCutoff stats <> title
                    (printf "%s: %d cells passed QC" (T.unpack $ input^.eid) (length stats'))
            return (cellQC, (input^.eid, stats'))
        savePlots output cellQCs 
            [ plotNumReads stats, plotTE stats, plotDupRate stats
            , plotMitoRate stats ]
        B.writeFile outputStat $ B.unlines $ map showStat $
            flip concatMap stats $ \(i, stat) -> map (\x ->
                x{_barcode = B.pack (T.unpack i) <> "+" <> _barcode x}) stat
        return outputStat

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
rmDup [] = return []
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
rmChrM _ [] = return ([], [])
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
        tss | geneStrand = geneLeft : map transLeft geneTranscripts
            | otherwise = geneRight : map transRight geneTranscripts

removeDoublet :: SCATACSeqConfig config
              => SCATACSeq S ( File '[NameSorted, Gzip] 'Bed
                             , File '[] 'Tsv )
              -> ReaderT config IO (SCATACSeq S (File '[NameSorted, Gzip] 'Bed, Int))
removeDoublet input = do
    dir <- asks _scatacseq_output_dir >>= getPath . (<> (asDir "/Bed"))
    teCutoff <- asks _scatacseq_te_cutoff
    let output = printf "%s/%s_rep%d_srt_filt_no_dblet.bed.gz" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . (\(bedFl, statFl) -> do
        stats <- readStats $ statFl^.location
        let cells = S.fromList $ map _barcode $ filter (passedQC teCutoff) stats
            f :: [BED] -> Bool
            f = (`S.member` cells) . fromJust . (^.name) . head
            sink = zipSinks lengthC $ concatC .| sinkFileBedGzip output
        (nCells, _) <- runResourceT $ runConduit $ streamBedGzip (bedFl^.location) .|
            groupBy ((==) `on` (^.name)) .| filterC f .| sink
        return (location .~ output $ emptyFile, nCells) )

detectDoublet :: SCATACSeqConfig config
              => SCATACSeq S (File tags 'Other, (a, b, File '[] 'Tsv ))
              -> ReaderT config IO (SCATACSeq S (File '[] 'Tsv))
detectDoublet input = do
    dir <- qcDir
    let output = dir <> "doublet_" <> T.unpack (input^.eid) <> "_rep" <>
            show (input^.replicates._1) <> ".tsv"
        outputPlot = dir <> "doublet_" <> T.unpack (input^.eid) <> "_rep" <>
            show (input^.replicates._1) <> ".html"
    input & replicates.traverse.files %%~ liftIO . (\(fl, (_,_,stat)) -> withTemp Nothing $ \tmp -> do
        shelly $ run_ "taiji-utils" ["doublet", T.pack $ fl^.location, T.pack tmp]
        [probs, threshold, sc, sim_sc] <- B.lines <$> B.readFile tmp
        let th = readDouble threshold
            dProbs = map readDouble $ B.words probs
            ds = map readDouble $ B.words sc
            ds_sim = map readDouble $ B.words sim_sc
            rate = fromIntegral (length $ filter (>0.5) dProbs) /
                fromIntegral (length dProbs) * 100 :: Double
        savePlots outputPlot [ mkHist ds th <> title (printf "doublet percentage: %.1f%%" rate)
            , mkHist ds_sim th ] []

        mat <- mkSpMatrix id $ fl^.location
        bcs <- runResourceT $ runConduit $ streamRows mat .| mapC fst .| sinkList
        let doubletScore = M.fromList $ zip bcs dProbs
        stats <- readStats $ stat^.location
        B.writeFile output $ B.unlines $ flip map stats $ \s ->
            let v = M.findWithDefault 0 (_barcode s) doubletScore
            in showStat s{_doublet_score = v}
        return $ location .~ output $ emptyFile
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

plotDupRate :: [(T.Text, [Stat])] -> E.EChart
plotDupRate stats = E.addAttr [jmacroE| {
    yAxis: {
        name: "duplication rate",
        nameLocation: "middle",
        nameGap: 50
    }
    }|] plt
  where
    plt = E.boxplot $ flip map stats $ \(nm, stat) -> 
        (nm, map _dup_rate stat)

plotMitoRate :: [(T.Text, [Stat])] -> E.EChart
plotMitoRate stats = E.boxplot $ flip map stats $ \(nm, stat) -> 
    (nm, map _mito_rate stat)

plotNumReads :: [(T.Text, [Stat])] -> E.EChart
plotNumReads stats = E.addAttr [jmacroE| {
    yAxis: {
        name: "log10(number of fragments)",
        nameLocation: "middle",
        nameGap: 50
    }
    }|] plt
  where
    plt = E.boxplot $ flip map stats $ \(nm, stat) -> 
        (nm, map (logBase 10 . fromIntegral . _uniq_reads) stat)

plotTE :: [(T.Text, [Stat])] -> E.EChart
plotTE stats = E.addAttr [jmacroE| {
    yAxis: {
        name: "TSS enrichment",
        nameLocation: "middle",
        nameGap: 50
    }
    }|] plt
  where
    plt = E.boxplot $ flip map stats $ \(nm, stat) -> 
        (nm, map _te stat)

plotCells :: Double -> [Stat] -> Vega
plotCells teCutoff input = plt <> axes <> vline <> hline <> scales
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
                    y: {value: `teCutoff`, scale: "y"},
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

--------------------------------------------------------------------------------
-- Cluster level QC
--------------------------------------------------------------------------------

{-
plotClusterQC :: SCATACSeqConfig config
              => ( FilePath
                 , SCATACSeq S
                    ( File '[] 'Tsv
                    , File '[] 'Other
                    , File '[Gzip] 'Other
                    , [(T.Text, File '[] 'Tsv)] ) )
              -> ReaderT config IO ()
plotClusterQC (idxFl, input) = do
    dir <- qcDir
    let output = dir <> T.unpack (input^.eid) <> "_cluster_qc.html"
    (geneList, plt1, plt2) <- plotDiffGene diff
    liftIO $ do
        stats <- M.fromList . map (\x -> (_barcode x, x)) <$>
            readStats (statFl^.location)
        cls <- decodeFile $ clFl^.location
        geneExpr <- M.fromList <$>
            readGeneExpr (map (B.pack . T.unpack) geneList) idxFl matFl
        let points = flip map cls $ \(CellCluster nm cells) ->
                (B.unpack nm, map _cell_3d cells)
            statViz = [ ("read depth", map (logBase 10 . fromIntegral . _uniq_reads) statOrdered)
                , ("TSS enrichment", map _te statOrdered)
                , ("duplication rate", map _dup_rate statOrdered)
                , ("mito rate", map _mito_rate statOrdered)
                , ("doublet score", map _doublet_score statOrdered) ]
            geneViz = zipWith (\x y -> (x, map fromIntegral y))
                (map T.unpack geneList) $ transpose $ concatMap (map (\x ->
                    M.findWithDefault undefined (_cell_barcode x) geneExpr) .
                    _cluster_member) cls
            statOrdered = concatMap (map
                (\x -> M.findWithDefault undefined (_cell_barcode x) stats) .
                _cluster_member) cls
        savePlots output []
            [ E.scatter3D' points (statViz <> geneViz) <> E.toolbox
            , plt1, plt2 ]
  where
    (statFl, clFl, matFl, diff) = input^.replicates._2.files

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

plotDiffGene :: SCATACSeqConfig config
             => [(T.Text, File '[] 'Tsv)]
             -> ReaderT config IO ([T.Text], E.EChart, E.EChart)
plotDiffGene inputs = do
    (cls, fdrs) <- liftIO $ fmap unzip $ forM inputs $ \(cl, fl) -> do
        fdr <- readFDR $ fl^.location
        return (cl, fdr)
    markers <- asks _scatacseq_marker_gene_list >>= \case
        Nothing -> return []
        Just fl -> liftIO $ readMarkers fl
    let genes = nubSort $ concatMap M.keys fdrs
        df1 = DF.mkDataFrame cls (map fst markers) $ flip map fdrs $ \fdr ->
            U.toList $ scale' $ U.fromList $ map snd $ enrichment markers fdr
        df2 = DF.mkDataFrame cls genes $ flip map fdrs $ \fdr ->
            map (\g -> M.findWithDefault 0 g fdr) genes
    return ( filter (`elem` concatMap snd markers) genes
           , mkHeatmap df1, mkHeatmap df2 )
  where
    mkHeatmap df = p <> E.toolbox <> E.option
        [jmacroE| { visualMap: { inRange: {color: `viridis`} } } |]
      where
        p = E.heatmap $ DF.orderDataFrame id df
    scale' xs | U.all (==0) xs = xs
              | otherwise = scale xs
    enrichment :: [(T.Text, [T.Text])] -> M.Map T.Text Double
               -> [(T.Text, Double)]
    enrichment markers fdr = map (second f) markers
      where
        f xs = foldl' (+) 0 (map (\k -> M.findWithDefault 0 k fdr) xs) /
            fromIntegral (length xs)
    readFDR :: FilePath -> IO (M.Map T.Text Double)
    readFDR fl = M.fromList . map snd . filter ((>2) . fst) .
        map (f . T.splitOn "\t") . T.lines <$> T.readFile fl
      where
        f (a:fld:_:q:_) = (read $ T.unpack fld :: Double, (a, read $ T.unpack q))
        f _ = undefined
    readMarkers :: FilePath -> IO [(T.Text, [T.Text])]
    readMarkers fl = M.toList . M.fromListWith (++) . map (f . T.splitOn "\t") .
        T.lines <$> T.readFile fl
      where
        f (a:b:_) = (b, [T.toUpper a])
        f _ = undefined
-}