{-# LANGUAGE DataKinds #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE QuasiQuotes       #-}
{-# LANGUAGE TemplateHaskell   #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE ViewPatterns #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Clustering.SnapTools
   ( snapPre
   , performSnap
   ) where

import qualified Data.ByteString.Char8 as B
import Bio.Seq.IO
import Shelly hiding (FilePath)
import qualified Data.Text as T
import System.IO.Temp (withTempFile)
import System.IO

import qualified Language.R                        as R
import           Language.R.QQ
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

performSnap :: SCATACSeqConfig config 
            => SCATACSeq S (File '[] 'Other)
            -> ReaderT config IO (SCATACSeq S (File '[] 'Tsv, File '[Gzip] 'Tsv))
performSnap input = input & replicates.traverse.files %%~ ( \fl -> do
    dir <- asks ((<> "/Snap/") . _scatacseq_output_dir) >>= getPath
    let output1 = printf "%s/%s_rep%d_snap_rownames.tsv" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
        output2 = printf "%s/%s_rep%d_snap.tsv.gz" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    liftIO $ snap output1 output2 (fl^.location) 
    return ( location .~ output1 $ emptyFile
           , location .~ output2 $ emptyFile )
    )

snap :: FilePath   -- ^ Row names
     -> FilePath   -- ^ Matrix
     -> FilePath   -- ^ Input
     -> IO ()
snap rownames mat input = R.runRegion $ do
    _ <- [r| library("SnapATAC")
        x.sp <- createSnap(file=input_hs, sample="SNAP")
        x.sp <- addBmatToSnap(x.sp, bin.size=5000, num.cores=1)
        x.sp <- makeBinary(x.sp, mat="bmat")
        x.sp <- runJDA( obj=x.sp, input.mat="bmat",
            bin.cov.zscore.lower=-2,
            bin.cov.zscore.upper=2,
            pc.num=30,
            norm.method="normOVE",
            max.var=5000,
            do.par=TRUE,
            ncell.chunk=1000,
            num.cores=10,
            seed.use=10,
            tmp.folder=tempdir()
        )
        write.table(cbind(x.sp@barcode, rowSums(x.sp@bmat)),
            file=rownames_hs, sep="\t", row.names=F, col.names=F, quote=F)
        gz1 <- gzfile(mat_hs, "w")
        write.table(x.sp@smat@dmat, gz1, sep="\t", row.names=F, col.names=F, quote=F)
        close(gz1)
    |]
    return ()