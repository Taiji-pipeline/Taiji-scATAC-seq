{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE QuasiQuotes #-}
{-# LANGUAGE DataKinds #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.QC
    ( plotStat
    , readStats
    , writeStats
    , showStat
    , decodeStat
    , streamTN5Insertion
    , rmAbnoramlFragment
    , rmChrM
    , rmDup
    , tssEnrichment
    , readTSS
    , detectDoublet
    , plotNumReads
    , plotDupRate
    , plotMitoRate
    , plotTE
    , plotDoubletScore
    ) where

import           Bio.Data.Bam
import           Bio.HTS
import Control.Monad.State.Strict
import Language.Javascript.JMacro
import Bio.Utils.Functions (slideAverage)
import           Bio.Data.Bed
import Bio.RealWorld.GENCODE
import qualified Data.ByteString.Char8 as B
import qualified Data.Map.Strict as M
import qualified Data.HashSet as S
import Data.Conduit.List (groupBy)
import qualified Data.Text as T
import Data.List.Ordered (nubSort)
import qualified Data.IntervalMap.Strict      as IM
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Unboxed.Mutable as UM
import Shelly (shelly, run_)

import Taiji.Prelude hiding (groupBy)
import Taiji.Pipeline.SC.ATACSeq.Types
import Taiji.Pipeline.SC.ATACSeq.Functions.Utils (extractBarcode)
import Taiji.Utils.Plot
import Taiji.Utils.Plot.Vega
import qualified Taiji.Utils.Plot.ECharts as E
import Taiji.Utils (mkSpMatrix, streamRows, readGenesValidated)

plotStat :: SCATACSeqConfig config
         => [SCATACSeq S (TagsAligned)]
         -> ReaderT config IO FilePath
plotStat [] = return ""
plotStat inputs = do
    dir <- qcDir
    passedQC <- getQCFunction
    teCutoff <- asks _scatacseq_te_cutoff 
    let output = dir <> "/QC.html"
        outputStat = dir <> "/QC_merged.tsv"
    liftIO $ do
        (cellQCs, stats) <- fmap unzip $ forM inputs $ \input -> do
            stats <- readStats $ input^.replicates._2.files._2.location
            let stats' = filter passedQC stats
                cellQC = plotCells teCutoff stats <> title
                    (printf "%s: %d cells passed QC" (T.unpack nm) (length stats'))
                nm = (input^.eid) <> "_" <> T.pack (show $ input^.replicates._1)
            return (cellQC, (nm, stats'))
        savePlots output cellQCs 
            [ plotNumReads stats, plotTE stats, plotDupRate stats
            , plotMitoRate stats ]
        writeStats outputStat $ flip concatMap stats $ \(i, stat) -> map (\x ->
            x{_barcode = B.pack (T.unpack i) <> "+" <> _barcode x}) stat
        return outputStat

readStats :: FilePath -> IO [Stat]
readStats = fmap (map decodeStat . tail . B.lines) . B.readFile

writeStats :: FilePath -> [Stat] -> IO ()
writeStats output stats = B.writeFile output $ B.unlines $ header : map showStat stats
  where
    header = "barcode\tduplication_rate\tchrM_rate\tTSSe\tunique_fragment\tdoublet_likelihood"

showStat :: Stat -> B.ByteString
showStat Stat{..} = B.intercalate "\t" $
    [ _barcode
    , maybe "NA" (B.pack . show) _dup_rate
    , maybe "NA" (B.pack . show) _mito_rate
    , B.pack $ show _te
    , B.pack $ show _uniq_reads
    , maybe "NA" (B.pack . show) _doublet_score ]
{-# INLINE showStat #-}

decodeStat :: B.ByteString -> Stat
decodeStat x = Stat bc
    (if dupRate == "NA" then Nothing else Just $ readDouble dupRate)
    (if mitoRate == "NA" then Nothing else Just $ readDouble mitoRate)
    (readDouble te) (readInt uniq) 
    (if ds == "NA" then Nothing else Just $ readDouble ds)
  where
    [bc, dupRate, mitoRate, te, uniq, ds] = B.split '\t' x
{-# INLINE decodeStat #-}

{-
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
-}

rmAbnoramlFragment:: [BAM] -> [BAM]
rmAbnoramlFragment = filter $ \bam ->
    let s = abs $ tLen bam in s > 30 && s < 1000
{-# INLINE rmAbnoramlFragment #-}

-- | Remove duplicates for reads originated from a single cell and convert
-- BAM to BED with TN5 correction.
-- The BED interval of the fragment is obtained by adjusting the BAM alignment
-- interval of the sequenced read-pair. The start of the interval is moved
-- forward by 4bp from a left-most alignment position and backward 5bp from the
-- right-most alignment position. The transposase cuts the two DNA strands with
-- a 9bp overhang, and adjusted positions represent the center point between
-- these cuts; this position is recorded as a cut site that represents a
-- chromatin accessibility event.
rmDup :: Bool -> BAMHeader -> [BAM] -> [BED]
rmDup pair hdr input' = map (\(_, bed, c) -> score .~ Just c $ bed) $ M.elems $
    M.fromListWith collapse $ map (\b -> (getKey b, (getSc b, toBed b, 1))) input
  where
    input = filter (\b ->
        let f = flag b
        in not pair || (isFirstSegment f && hasMultiSegments f && not (isNextUnmapped f))) input'
    collapse (sc1, b1, c1) (sc2, b2, c2) | sc1 > sc2 = (sc1, b1, c1+c2)
                                         | otherwise = (sc2, b2, c1+c2)
    getKey b | pair = pairKey
             | otherwise = singleKey
      where
        (singleKey, pairKey) = makeKey (const Nothing) b
    getSc b = case queryAuxData ('m', 's') b of
        Just (AuxInt x) -> x + fromJust (sumQual 15 b)
        _ -> fromJust $ sumQual 15 b
    toBed b | pair = adjust $ fromJust $ bamToFragment hdr b
            | otherwise = adjust $ fromJust $ bamToBed hdr b
      where
        adjust x = name %~ fmap extractBarcode $ chromStart %~ (+4) $ chromEnd %~ (subtract 5) $ x
{-# INLINE rmDup #-}

-- | Remove chrM reads.
rmChrM :: [BED] -> ([BED], [BED], Double)
rmChrM [] = ([], [], 0)
rmChrM input = (output, mito, chrMRate)
  where
    (output, mito) = partition notChrM input
    chrMRate = fromIntegral (length mito) / fromIntegral (length input)
    notChrM x = let chr = x^.chrom in chr /= "chrM" && chr /= "M"
{-# INLINE rmChrM #-}

tssEnrichment :: BEDTree (Int, Bool)   -- ^ TSS
              -> [BED]   -- ^ Valid input can be fragments or single-end reads
              -> Double
tssEnrichment regions beds = U.maximum $ normalize $ U.create $ do
            v <- UM.replicate 4000 0
            forM_ (concatMap getCutSite beds) $ \cutsite ->
                forM_ (IM.elems $ intersecting regions cutsite) $ \(x, str) ->
                    UM.unsafeModify v (+1) $ if str
                        then cutsite^.chromStart - (x - 2000)
                        else 3999 - (x + 2000 - cutsite^.chromStart)
            return v
  where
    normalize vec = slideAverage 5 $ U.map (/bk) vec
      where
        bk = (U.sum (U.take 100 vec) + U.sum (U.drop 3900 vec)) / 200 + 0.1
    getCutSite bed = case bed^.strand of
        -- This is a fragment
        Nothing -> [left, right]
        -- Single forward strand
        Just True -> [left]
        -- Single reverse strand
        Just False -> [right]
      where
        left = let i = bed^.chromStart + 75 in BED3 (bed^.chrom) i $ i+1
        right = let i = bed^.chromEnd - 76 in BED3 (bed^.chrom) i $ i+1
{-# INLINE tssEnrichment #-}

readTSS :: FilePath -> IO (BEDTree (Int, Bool))
readTSS = fmap (bedToTree const . concatMap fn) . readGenesValidated
  where
    fn Gene{..} = map g $ nubSort tss
      where
        g x = (BED3 geneChrom (x - 2000) (x + 2000), (x, geneStrand))
        tss | geneStrand = geneLeft : map transLeft geneTranscripts
            | otherwise = geneRight : map transRight geneTranscripts
{-# INLINE readTSS #-}

-- | Filter QC-failed cells and convert fragments into TN5 insertions.
streamTN5Insertion :: ( unpair ~ '[NameSorted, Gzip]
                      , paired ~ '[NameSorted, PairedEnd, Gzip] )
                   => (Stat -> Bool)
                   -> (Either (File unpair 'Bed) (File paired 'Bed), File '[] 'Tsv)
                   -> ConduitT () [BED] (ResourceT IO) ()
streamTN5Insertion passedQC (input, qcFl) = do
    stats <- liftIO $ readStats $ qcFl^.location
    let cells = S.fromList $ map _barcode $ filter passedQC stats
    either (f cells . SomeFile) (f cells . SomeFile) input
  where
    f cells (SomeFile fl) = streamBedGzip (fl^.location) .|
        groupBy ((==) `on` (^.name)) .|
        filterC ((`S.member` cells) . fromJust . (^.name) . head) .|
        mapC (concatMap getCutSite . filter notChrM)
      where
        getCutSite x | fl `hasTag` PairedEnd = [left, right]
                     | otherwise = [x]
          where
            left = let i = x^.chromStart in BED (x^.chrom) i (i+1) (x^.name) Nothing (Just True)
            right = let i = x^.chromEnd - 1 in BED (x^.chrom) i (i+1) (x^.name) Nothing (Just False)
    notChrM x = let chr = x^.chrom in chr /= "chrM" && chr /= "M"
{-# INLINE streamTN5Insertion #-}

detectDoublet :: SCATACSeqConfig config
              => SCATACSeq S (File tags 'Matrix, TagsAligned)
              -> ReaderT config IO (SCATACSeq S TagsAligned)
detectDoublet input = do
    dir <- qcDir
    let outputPlot = dir <> "doublet_" <> T.unpack (input^.eid) <> "_rep" <>
            show (input^.replicates._1) <> ".html"
    input & replicates.traverse.files %%~ liftIO . (\(fl, (tags,stat)) -> withTemp Nothing $ \tmp -> do
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
        writeStats (stat^.location) $ flip map stats $ \s ->
            let v = M.findWithDefault 0 (_barcode s) doubletScore
            in s{_doublet_score = Just v}
        return (tags, stat)
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
        (nm, map (fromMaybe 0 . _dup_rate) stat)

plotMitoRate :: [(T.Text, [Stat])] -> E.EChart
plotMitoRate stats = E.boxplot $ flip map stats $ \(nm, stat) -> 
    (nm, map (fromMaybe 0 . _mito_rate) stat)

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

plotDoubletScore :: [(T.Text, [Stat])] -> E.EChart
plotDoubletScore stats = E.addAttr [jmacroE| {
    yAxis: {
        name: "probability of doublet",
        nameLocation: "middle",
        nameGap: 50
    }
    }|] plt
  where
    plt = E.boxplot $ flip map stats $ \(nm, stat) -> 
        (nm, map (fromMaybe 0 . _doublet_score) stat)

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
