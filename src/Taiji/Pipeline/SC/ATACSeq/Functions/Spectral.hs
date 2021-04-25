{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE GADTs #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Spectral
    ( getSpectral
    , nystromExtend
    , mergeResults
    , chunksInput
    ) where

import Data.Binary (decodeFile)
import Data.Conduit.Zlib (multiple, ungzip, gzip)
import qualified Data.Text as T
import System.Random.MWC
import Shelly
import Data.Hashable (hash)

import Taiji.Prelude
import Taiji.Utils
import Taiji.Pipeline.SC.ATACSeq.Types

chunksInput :: [File '[Gzip] 'Matrix]
            -> IO [([File '[Gzip] 'Matrix], Int)]
chunksInput mats = do
    sizes <- mapM (fmap _num_row . mkSpMatrix id . (^.location)) mats
    return $ go 0 $ zip3 (scanl1 (+) sizes) sizes mats
  where
    chunkSize = 50000
    go _ [] = []
    go start xs = let (res, (next, rest)) = takeN start xs in res : go next rest
    takeN start xs = ( (map (^._3) res, start)
                     , (nextStart, map f $ if nextStart == 0 then rest else last res : rest) )
      where
        f (a,b,c) = (a - chunkSize, b, c)
        nextStart = let (acc, n, _) = last res in if acc >= chunkSize then 0 else n - (acc - chunkSize)
        (res, rest) = partition ((<=chunkSize) . (^._1)) xs

getSpectral :: SCATACSeqConfig config
            => Int     -- ^ Number of samples
            -> ( [ SCATACSeq S ( File '[RowName, Gzip] 'Tsv
                               , File '[ColumnName, Gzip] 'Tsv
                               , File '[Gzip] 'Matrix ) ]
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
              => ( ([File '[Gzip] 'Matrix], Int)
                 , FilePath
                 , FilePath )
              -> ReaderT config IO (File '[Gzip] 'Tsv)
nystromExtend ((matFls, start), modelFl, cidxFl) = do
    tmpdir <- asks _scatacseq_tmp_dir
    dir <- asks ((<> "/temp/Nystrom/") . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%d.txt.gz" dir (hash (map (^.location) matFls, start))
    liftIO $ withTemp tmpdir $ \tmp -> do
        cidx <- decodeFile cidxFl
        ( deleteCols cidx . takeRows chunkSize . dropRows start .
            concatMatrix <$> mapM (mkSpMatrix id . (^.location)) matFls ) >>=
                saveMatrix tmp id
        shelly $ run_ "taiji-utils" $ ["predict", T.pack output,
            "--model", T.pack modelFl, "--input", T.pack tmp]
        return $ location .~ output $ emptyFile
  where
    chunkSize = 50000

mergeResults :: SCATACSeqConfig config
             => ( Maybe (Either FilePath FilePath)
                , [ SCATACSeq S ( File '[RowName, Gzip] 'Tsv
                                , File '[ColumnName, Gzip] 'Tsv
                                , File '[Gzip] 'Matrix ) ]
                , [File '[Gzip] 'Tsv] )
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
            mergeFiles (map (^.location) fls2) .| gzip .| sinkFile output
    return $ Just $ head fls1 & eid .~ "Merged" & replicates._2.files .~
        (location .~ rownames $ emptyFile, location .~ output $ emptyFile)
  where
    mergeFiles inFls = mapM_ source inFls .| unlinesAsciiC
      where
        source x = sourceFile x .| multiple ungzip .| linesUnboundedAsciiC