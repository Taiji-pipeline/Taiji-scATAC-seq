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
import Bio.Pipeline
import Shelly hiding (FilePath)
import qualified Data.Text as T
import System.IO.Temp (withTempFile)
import System.IO
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
                , "--bin-size-list", "5000", "50000" ]
            return $ location .~ output $ emptyFile )

getClusters :: SCATACSeq S (File '[] 'Other)
            -> IO (SCATACSeq S [CellCluster])
getClusters input = input & replicates.traverse.files %%~
    (\fl -> snap (fl^.location))

snap :: FilePath   -- ^ BED File
     -> IO [CellCluster]
snap input = R.runRegion $ do
    _ <- [r| library("SnapATAC")
        x.sp <- createSnap(file=input_hs, sample="SNAP")
        x.sp <- addBmatToSnap(x.sp, bin.size=50000, num.cores=1)
        x.sp <- makeBinary(x.sp, mat="bmat")
        x.sp <- filterBins( x.sp, low.threshold=-1.5
            , high.threshold=1.5, mat="bmat" )
        x.sp <- runJaccard( obj = x.sp, tmp.folder=tempdir(),
            mat = "bmat", max.var=2000, ncell.chunk=1000,
            do.par=FALSE, num.cores=1, seed.use=10 )
        x.sp <- runNormJaccard( obj = x.sp, tmp.folder=tempdir(),
            ncell.chunk=1000, method="normOVE", row.center=TRUE,
            row.scale=TRUE, low.threshold=-5, high.threshold=5,
            num.cores=1, seed.use=10 )
        x.sp <- runDimReduct( x.sp, pc.num=50, input.mat="jmat",
            method="svd", center=TRUE, scale=FALSE, seed.use=10 )
        x.sp <- runKNN( obj=x.sp, pca.dims=1:30, weight.by.sd=T,
            k=15 );
        x.sp <- runCluster(obj=x.sp, tmp.folder=tempdir(),
            louvain.lib="leiden", seed.use=10, resolution=1)
        barcodes <<- as.character(x.sp@metaData[,1])
        member <<- x.sp@cluster

        # Viz
        x.sp <- runViz( obj=x.sp, tmp.folder=tempdir(), dims=2,
            pca.dims=1:30, weight.by.sd=T, method="umap", Y.init=NULL,
            seed.use=10, num.cores=1 )
        dataViz <- as.data.frame(obj@umap)
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
        flip map cells $ \(bc, _, x, y) -> Cell bc x y