{-# LANGUAGE DataKinds #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE QuasiQuotes       #-}
{-# LANGUAGE TemplateHaskell   #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE ViewPatterns #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Clustering.SnapTools
   ( snapPre
   , performSnap
   , mkSnapMat
   ) where

import qualified Data.ByteString.Char8 as B
import Bio.Seq.IO
import Bio.Utils.Misc (readInt)
import Data.ByteString.Lex.Integral (packDecimal)
import Bio.Data.Bed
import Shelly hiding (FilePath)
import qualified Data.Text as T
import System.IO.Temp (withTempFile)
import System.IO
import Data.Conduit.Internal (zipSources)

import qualified Language.R                        as R
import           Language.R.QQ
import Language.R.HExp

import Taiji.Prelude
import Taiji.Pipeline.SC.ATACSeq.Types
import Taiji.Pipeline.SC.ATACSeq.Functions.Utils

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




mkSnapMat :: SCATACSeqConfig config 
          => SCATACSeq S ((File tags 'Bed, File tags 'Bed, Int), File t 'Other)
          -> ReaderT config IO (SCATACSeq S (File '[] 'Other))
mkSnapMat input = input & replicates.traverse.files %%~ ( \((_,x,_), y) -> do
    dir <- asks ((<> "/Snap/") . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_snap.rds" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    liftIO $ snap' output (x^.location) (y^.location)
    return $ location .~ output $ emptyFile
    )


snap' :: FilePath -- ^ output
      -> FilePath   -- ^ column names
      -> FilePath   -- ^ Sparse Matrix
      -> IO ()
snap' output nameFl matFl = withTempDir (Just "./") $ \tmpdir -> do
    regions <- runResourceT $ runConduit $ streamBedGzip nameFl .| mapC f .| sinkList
    let ridx = tmpdir <> "/ridx.txt"
        cidx = tmpdir <> "/cidx.txt"
        vals = tmpdir <> "/vals.txt"
    mat <- mkSpMatrix readInt matFl
    rows <- parseSpMat ridx cidx vals mat
    R.runRegion $ do
        _ <- [r| library("Matrix")
                 i <- scan(ridx_hs, integer())
                 j <- scan(cidx_hs, integer())
                 vals <- scan(vals_hs, integer())
                 mat <- sparseMatrix(i,j,x=vals)
                 colnames(mat) <- regions_hs
                 rownames(mat) <- rows_hs
                 saveRDS(mat, file=output_hs)
        |]
        return ()
  where
    f :: BED3 -> String
    f bed = B.unpack (bed^.chrom) <> ":" <> show (bed^.chromStart) <> "-" <>
        show (bed^.chromEnd)

parseSpMat :: FilePath  -- ^ row index
           -> FilePath  -- ^ col index
           -> FilePath  -- ^ value
           -> SpMatrix Int
           -> IO [String]
parseSpMat ridx cidx vals mat = do
    (res,_,_,_) <- runResourceT $ runConduit $
        zipSources (iterateC succ 0) (streamRows mat) .| mapC f .|
        getZipSink ((,,,) <$> ZipSink sink1 <*> ZipSink sink2 <*> ZipSink sink3 <*> ZipSink sink4)
    return res
  where
    sink1 = mapC fst .| sinkList
    sink2 = concatMapC (map (fromJust . packDecimal) . (^._1) . snd) .| intersperseC " " .| sinkFile ridx
    sink3 = concatMapC (map (fromJust . packDecimal) . (^._2) . snd) .| intersperseC " " .| sinkFile cidx
    sink4 = concatMapC (map (fromJust . packDecimal) . (^._3) . snd) .| intersperseC " " .| sinkFile vals
    f (i, (bc, xs)) = (B.unpack bc
        , unzip3 $ map (\(i, (j,x)) -> (i,j,if x > 0 then 1 else 0)) $ zip (repeat i) xs)