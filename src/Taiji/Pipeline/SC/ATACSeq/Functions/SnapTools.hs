{-# LANGUAGE DataKinds #-}
{-# LANGUAGE OverloadedStrings #-}
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
                , "--bin-size-list 5000 50000" ]
            return $ location .~ output $ emptyFile )