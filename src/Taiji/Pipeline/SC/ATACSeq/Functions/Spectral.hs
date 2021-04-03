{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE GADTs #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Spectral
    ( getSpectral
    , nystromExtend
    , mergeResults
    ) where

import Data.Binary (decodeFile)
import Data.Conduit.Zlib (multiple, ungzip, gzip)
import qualified Data.Text as T
import System.Random.MWC
import Shelly

import Taiji.Prelude
import Taiji.Utils
import Taiji.Pipeline.SC.ATACSeq.Types

getSpectral :: SCATACSeqConfig config
            => Int     -- ^ Number of samples
            -> ( [ SCATACSeq S ( File '[RowName, Gzip] 'Tsv
                               , File '[ColumnName, Gzip] 'Tsv
                               , File '[Gzip] 'Other ) ]
               , Maybe FilePath )
            -> ReaderT config IO (Maybe (Either FilePath FilePath))
getSpectral sampleSize (input, Just cidxFl) = do
    dir <- asks ((<> "/Spectral/") . _scatacseq_output_dir) >>= getPath
    let output = dir <> "model.pbz2"
    tmpdir <- asks _scatacseq_tmp_dir
    liftIO $ withTemp tmpdir $ \tmp -> do
        cidx <- decodeFile cidxFl
        mat <- fmap (deleteCols cidx . concatMatrix) $ forM input $ \x ->
            mkSpMatrix id $ x^.replicates._2.files._3.location
        if _num_row mat <= sampleSize
            then do
                saveMatrix tmp id mat
                shelly $ run_ "taiji-utils" ["fit", T.pack tmp, T.pack output]
                return $ Just $ Left output
            else do
                create >>= sampleRows sampleSize mat >>= saveMatrix tmp id
                let rate = fromIntegral sampleSize / fromIntegral (_num_row mat) :: Double
                shelly $ run_ "taiji-utils" ["fit", T.pack tmp, T.pack output,
                    "--sampling-rate", T.pack $ show rate]
                return $ Just $ Right output
getSpectral _ _ = return Nothing

nystromExtend :: SCATACSeqConfig config
              => ( SCATACSeq S ( File '[RowName, Gzip] 'Tsv
                               , File '[ColumnName, Gzip] 'Tsv
                               , File '[Gzip] 'Other )
                 , FilePath
                 , FilePath )
              -> ReaderT config IO (SCATACSeq S (File '[Gzip] 'Tsv))
nystromExtend (input, modelFl, cidxFl) = do
    tmpdir <- asks _scatacseq_tmp_dir
    dir <- asks ((<> "/Spectral/Sample/") . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_nystrom.txt.gz" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . (\(_, _, matFl) -> withTemp tmpdir $ \tmp -> do
        cidx <- decodeFile cidxFl
        (deleteCols cidx <$> mkSpMatrix id (matFl^.location)) >>= saveMatrix tmp id
        shelly $ run_ "taiji-utils" $ ["predict", T.pack output,
            "--model", T.pack modelFl, "--input", T.pack tmp]
        return $ location .~ output $ emptyFile
        )

mergeResults :: SCATACSeqConfig config
             => ( Maybe (Either FilePath FilePath)
                , [ SCATACSeq S ( File '[RowName, Gzip] 'Tsv
                                , File '[ColumnName, Gzip] 'Tsv
                                , File '[Gzip] 'Other ) ]
                , [SCATACSeq S (File '[Gzip] 'Tsv)] )
             -> ReaderT config IO (Maybe (SCATACSeq S (File '[] 'Tsv, File '[Gzip] 'Tsv)))
mergeResults (Nothing, _, _) = return Nothing
mergeResults (Just modelFl, fls1, fls2) = do
    dir <- asks ((<> "/Spectral/") . _scatacseq_output_dir) >>= getPath
    let rownames = dir <> "Merged_rownames.txt"
        output = dir <> "Merged_reduced_dims.txt.gz"
    liftIO $ runResourceT $ runConduit $
        mergeFiles (map (^.replicates._2.files._1.location) fls1) .| sinkFile rownames
    liftIO $ case modelFl of
        Left model -> shelly $ run_ "taiji-utils" $ ["predict", T.pack output,
            "--model", T.pack model]
        Right _ -> runResourceT $ runConduit $
            mergeFiles (map (^.replicates._2.files.location) fls2) .| gzip .| sinkFile output
    return $ Just $ head fls1 & eid .~ "Merged" & replicates._2.files .~
        (location .~ rownames $ emptyFile, location .~ output $ emptyFile)
  where
    mergeFiles inFls = mapM_ source inFls .| unlinesAsciiC
      where
        source x = sourceFile x .| multiple ungzip .| linesUnboundedAsciiC