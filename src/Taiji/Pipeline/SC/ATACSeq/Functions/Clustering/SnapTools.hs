{-# LANGUAGE DataKinds #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE QuasiQuotes       #-}
{-# LANGUAGE TemplateHaskell   #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE ViewPatterns #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Clustering.SnapTools
   ( snapPre
   , getClusters
   ) where

import qualified Data.ByteString.Char8 as B
import qualified Data.ByteString as BS
import Bio.Seq.IO
import Shelly hiding (FilePath)
import qualified Data.Text as T
import System.IO.Temp (withTempFile)
import System.IO
import Data.Binary
import Data.Int (Int32)

import qualified Language.R                        as R
import           Language.R.QQ
import qualified Data.Vector.SEXP as V
import Language.R.HExp

import Taiji.Prelude
import Taiji.Pipeline.SC.ATACSeq.Types

snapPre :: SCATACSeqConfig config
        => SCATACSeq S (File '[NameSorted, Gzip] 'Bed)
        -> ReaderT config IO (SCATACSeq S (File '[] 'Other))
snapPre input = do
    dir <- asks ((<> "/Snap") . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d.snap" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    genome <- asks (fromJust . _scatacseq_genome_index)
    chrSizes <- liftIO $ withGenome genome $ return . getChrSizes
    input & replicates.traverse.files %%~ ( \fl -> liftIO $
        withTempFile "./" "tmp_chrsize" $ \tmpChr h -> do
            B.hPutStr h $ B.unlines $
                map (\(a,b) -> a <> "\t" <> B.pack (show b)) chrSizes
            hClose h
            shelly $ run_ "snaptools" [ "snap-pre"
                , "--input-file=" <> T.pack (fl^.location)
                , "--output-snap=" <> T.pack output
                , "--genome-name=hg19"
                , "--genome-size=" <> T.pack tmpChr
                , "--min-flen=0"
                , "--max-flen=1000"
                , "--overwrite=True"
                , "--min-cov=100" ]
            shelly $ run_ "snaptools" [ "snap-add-bmat"
                , "--snap-file=" <> T.pack output
                , "--bin-size-list", "5000" ]
            return $ location .~ output $ emptyFile )

getClusters :: SCATACSeqConfig config 
            => SCATACSeq S (File '[] 'Other)
            -> ReaderT config IO (SCATACSeq S (File '[] 'Other))
getClusters input = input & replicates.traverse.files %%~ ( \fl -> do
    dir <- asks ((<> "/Snap/") . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_clusters.bin" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    liftIO $ snap (fl^.location) >>= encodeFile output
    return $ location .~ output $ emptyFile )

snap :: FilePath   -- ^ BED File
     -> IO [CellCluster]
snap input = R.runRegion $ do
    _ <- [r| library("SnapATAC")
             library("leiden")
        x.sp <- createSnap(file=input_hs, sample="SNAP")
        x.sp <- addBmatToSnap(x.sp, bin.size=5000, num.cores=1)
        x.sp <- makeBinary(x.sp, mat="bmat")

        x.sp <- runJDA( obj=x.sp, input.mat="bmat",
            bin.cov.zscore.lower=-2,
            bin.cov.zscore.upper=2,
            pc.num=50,
            norm.method="normOVE",
            max.var=5000,
            do.par=TRUE,
            ncell.chunk=1000,
            num.cores=10,
            seed.use=10,
            tmp.folder=tempdir()
        )
        x.sp <- runKNN( obj=x.sp,
            pca.dims=2:40,
            weight.by.sd=TRUE,
            k=15
        )
        x.sp <- runCluster(obj=x.sp, tmp.folder=tempdir(),
            louvain.lib="R-igraph", seed.use=10, resolution=1)
        barcodes <<- as.character(x.sp@metaData[,1])
        member <<- x.sp@cluster

        # Viz
        x.sp <- runViz( obj=x.sp, tmp.folder=tempdir(), dims=2,
            pca.dims=2:40, weight.by.sd=T, method="Rtsne", Y.init=NULL,
            seed.use=10, num.cores=1 )
        dataViz <- as.data.frame(x.sp@tsne)
        vizX <<- dataViz[,1]
        vizY <<- dataViz[,2]
    |]
    membership <- [r| member |] >>= ( \x -> R.unSomeSEXP x $
        \(hexp -> Int v) -> return $ V.toList v )
    bc <- [r| barcodes |] >>= ( \x -> R.unSomeSEXP x $
        \(hexp -> String v) -> return $ map
        (\(hexp -> Char c) -> BS.pack $ V.toList c) $ V.toList v )
    coordX <- [r| vizX |] >>= ( \x -> R.unSomeSEXP x $
        \(hexp -> Real v) -> return $ V.toList v )
    coordY <- [r| vizY |] >>= ( \x -> R.unSomeSEXP x $
        \(hexp -> Real v) -> return $ V.toList v )
    return $ zipWith f [0..] $ groupBy ((==) `on` (^._2)) $
        sortBy (comparing (^._2)) $ zip4 bc (membership :: [Int32]) coordX coordY
  where
    f i cells = CellCluster (B.pack $ "C" ++ show (i :: Int)) $
        flip map cells $ \(bc, _, x, y) -> Cell 1 x y 0 bc 0