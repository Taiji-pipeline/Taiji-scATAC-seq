{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE DataKinds #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Clustering.DiffusionMap (performDM) where 

import Data.ByteString.Lex.Integral (packDecimal)
import Data.Singletons.Prelude (Elem)
import Bio.Utils.Misc (readInt)
import qualified Data.Text as T
import Shelly (shelly, run_)

import Taiji.Prelude
import Taiji.Pipeline.SC.ATACSeq.Types
import Taiji.Pipeline.SC.ATACSeq.Functions.Utils

performDM :: (Elem 'Gzip tags ~ 'True, SCATACSeqConfig config)
          => FilePath  -- ^ directory
          -> SCATACSeq S (File tags 'Other)
          -> ReaderT config IO (SCATACSeq S (File '[] 'Tsv, File '[Gzip] 'Tsv))
performDM prefix input = do
    dir <- asks ((<> asDir prefix) . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_dm.tsv.gz" dir
            (T.unpack $ input^.eid) (input^.replicates._1)
        rownames = printf "%s/%s_rep%d_dm.rownames.txt" dir
            (T.unpack $ input^.eid) (input^.replicates._1)
    input & replicates.traversed.files %%~ liftIO . ( \fl -> do
        diffusionMap output fl
        sp <- mkSpMatrix readInt $ fl^.location
        runResourceT $ runConduit $
            streamRows sp .| mapC f .| unlinesAsciiC .| sinkFile rownames
        return ( location .~ rownames $ emptyFile
               , location .~ output $ emptyFile ) )
  where
    f (nm, xs) = nm <> "\t" <> fromJust (packDecimal $ foldl1' (+) $ map snd xs)


diffusionMap :: Elem 'Gzip tags ~ 'True
             => FilePath
             -> File tags 'Other
             -> IO ()
diffusionMap output input = shelly $ run_ "sc_utils"
    ["reduce", "--method", "dm", T.pack $ input^.location, T.pack output]