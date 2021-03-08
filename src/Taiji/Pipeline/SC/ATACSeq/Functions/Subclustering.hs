{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE DataKinds #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Subclustering
    ( subSpectral ) where

import qualified Data.ByteString.Char8 as B
import Data.Binary (encodeFile, decodeFile)
import Bio.Utils.Functions (scale)
import Data.Conduit.Zlib (multiple, ungzip, gzip)
import qualified Data.Conduit.List as CL
import qualified Data.Text as T
import qualified Data.HashSet as S
import qualified Data.HashMap.Strict as M
import Data.Singletons.Prelude (Elem)
import Bio.Data.Bed
import Control.Arrow (first, (&&&))
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed.Mutable as UM
import qualified Data.Vector.Unboxed as U
import System.IO
import Data.List.Ordered (nubSort)
import Shelly (shelly, run_, escaping)
import Control.Workflow
   
import Taiji.Pipeline.SC.ATACSeq.Functions.QC
import Taiji.Prelude
import Taiji.Pipeline.SC.ATACSeq.Types
import qualified Taiji.Utils.DataFrame as DF
import Taiji.Utils.Plot
import Taiji.Utils
import Taiji.Utils.Plot.ECharts

-- | Reduce dimensionality using spectral clustering
subSpectral :: SCATACSeqConfig config
            => SCATACSeq S (File tags 'Other)
            -> ReaderT config IO (SCATACSeq S (File '[] 'Tsv, File '[Gzip] 'Tsv))
subSpectral input = do
    dir <- asks ((<> "/Subcluster/Spectral/") . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_spectral.tsv.gz" dir (T.unpack $ input^.eid)
        outputIdx = printf "%s/%s_rownames.tsv" dir (T.unpack $ input^.eid)
    input & replicates.traversed.files %%~ ( \fl -> do
        asks _scatacseq_batch_info >>= \case
            Nothing -> liftIO $ shelly $ run_ "taiji-utils" $ ["reduce", T.pack $ fl^.location,
                T.pack output, "--seed", "23948"]
            Just batchFl -> liftIO $ do
                idToBatchMap <- M.fromListWith undefined <$> readBatchInfo batchFl
                let f x = let (i, r) = B.breakEnd (=='_') x
                        in case M.lookup (B.init i) idToBatchMap of
                                Nothing -> Nothing
                                Just (l, g) -> Just (l <> r, g)
                mat <- mkSpMatrix id $ fl^.location
                barcodes <- runResourceT $ runConduit $ streamRows mat .| mapC fst .| sinkList
                B.writeFile outputIdx $ B.unlines barcodes
                let labels = map (f . B.init . fst . B.breakEnd (=='+')) barcodes
                if (all isNothing labels)
                    then shelly $ run_ "taiji-utils" $ ["reduce", T.pack $ fl^.location,
                        T.pack output, "--seed", "23948"]
                    else withTempDir (Just "./") $ \tmpdir -> do
                        let tmp = tmpdir <> "/tmp.tsv.gz"
                        shelly $ run_ "taiji-utils" $ ["reduce", T.pack $ fl^.location,
                            T.pack tmp, "--seed", "23948"]
                        readData tmp >>= batchCorrect labels >>= writeData output
        return ( location .~ outputIdx $ emptyFile
               , location .~ output $ emptyFile )
        )
  where
    readData fl = runResourceT $ runConduit $
        sourceFile fl .| multiple ungzip .| linesUnboundedAsciiC .|
        mapC (U.fromList . map readDouble . B.split '\t') .| sinkVector
    writeData output vec = runResourceT $ runConduit $ yieldMany vec .|
        mapC (B.intercalate "\t" . map toShortest . U.toList) .|
        unlinesAsciiC .| gzip .| sinkFile output