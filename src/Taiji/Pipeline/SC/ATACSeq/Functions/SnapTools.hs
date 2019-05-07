{-# LANGUAGE DataKinds #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE QuasiQuotes       #-}
{-# LANGUAGE TemplateHaskell   #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.SnapTools where

import qualified Data.ByteString.Char8 as B
import Bio.Seq.IO
import Bio.Pipeline
import Data.Maybe
import Scientific.Workflow
import Control.Monad.Reader (asks, liftIO)
import Text.Printf (printf)
import Bio.Data.Experiment
import Shelly hiding (FilePath)
import qualified Data.Text as T
import System.IO.Temp (withTempFile)
import System.IO
import Control.Lens

import qualified Language.R                        as R
import           Language.R.QQ

import Taiji.Pipeline.SC.ATACSeq.Types

snapPre :: SCATACSeqConfig config
        => SCATACSeq S (File '[NameSorted, Gzip] 'Bed)
        -> WorkflowConfig config (SCATACSeq S (File '[] 'Other))
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

snap :: FilePath -> IO ()
snap input = R.runRegion $ do
    [r| library(SnapATAC)
        x.sp <- createSnap(file=input_hs);
        x.sp <- addBmatToSnap(x.sp, bin.size=50000, num.cores=1);
        x.sp <- makeBinary(x.sp, mat="bmat");
        x.sp <- filterBins( x.sp, low.threshold=-1.5
            , high.threshold=1.5, mat="bmat" );
        x.sp <- runJaccard( obj = x.sp, tmp.folder=tempdir(),
            mat = "bmat", max.var=2000, ncell.chunk=1000,
            do.par=FALSE, num.cores=1, seed.use=10 );
        x.sp <- runNormJaccard( obj = x.sp, tmp.folder=tempdir(),
            ncell.chunk=1000, method="normOVE", row.center=TRUE,
            row.scale=TRUE, low.threshold=-5, high.threshold=5,
            num.cores=1, seed.use=10 );
        x.sp <- runDimReduct( x.sp, pc.num=50, input.mat="jmat",
            method="svd", center=TRUE, scale=FALSE, seed.use=10 );
        x.sp <- runKNN( obj=x.sp, pca.dims=1:30, weight.by.sd=FALSE,
            k=15 );
        x.sp <- runCluster( obj=x.sp, tmp.folder=tempdir(),
            louvain.lib="R-igraph", path.to.snaptools=NULL, seed.use=10 );
    |]
    return ()