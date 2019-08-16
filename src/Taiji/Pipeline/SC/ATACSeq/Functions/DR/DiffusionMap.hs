{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE DataKinds #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.DR.DiffusionMap (performDM) where 

import Data.ByteString.Lex.Integral (packDecimal)
import Data.Singletons.Prelude (Elem)
import Bio.Utils.Functions (scale)
import qualified Data.Vector.Unboxed.Mutable as UM
import qualified Data.Vector.Unboxed as U
import Bio.Utils.Misc (readInt)
import qualified Data.Text as T
import Shelly (shelly, run_)

import Taiji.Prelude
import Taiji.Pipeline.SC.ATACSeq.Types
import Taiji.Pipeline.SC.ATACSeq.Functions.Utils

performDM :: (Elem 'Gzip tags ~ 'True, SCATACSeqConfig config)
          => FilePath  -- ^ directory
          -> SCATACSeq S (File tags 'Other)
          -> ReaderT config IO (SCATACSeq S (File '[] 'Tsv, [File '[Gzip] 'Tsv]))
performDM prefix input = do
    dir <- asks ((<> asDir prefix) . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_dm.tsv.gz" dir
            (T.unpack $ input^.eid) (input^.replicates._1)
        rownames = printf "%s/%s_rep%d_dm.rownames.txt" dir
            (T.unpack $ input^.eid) (input^.replicates._1)
    input & replicates.traversed.files %%~ liftIO . ( \fl -> withTemp Nothing $ \tmp -> do
        sp <- mkSpMatrix readInt $ fl^.location

        -- filtering
        {-
        vec <- UM.replicate (_num_col sp) 0
        runResourceT $ runConduit $ streamRows sp .| concatMapC snd .|
            mapC fst .| mapM_C (UM.unsafeModify vec (+1))
        v <- scale <$> U.unsafeFreeze vec
        let idx = U.toList $ U.imapMaybe
                (\i x -> if x < -2 || x > 2 then Just i else Nothing) v
        filterCols tmp idx $ fl^.location

        diffusionMap output tmp
        -}
        diffusionMap output $ fl^.location

        runResourceT $ runConduit $
            streamRows sp .| mapC f .| unlinesAsciiC .| sinkFile rownames

        let o = [ location .~ output $ emptyFile ]
        return ( location .~ rownames $ emptyFile
               , o ) )
  where
    f (nm, xs) = nm <> "\t" <> fromJust (packDecimal $ foldl1' (+) $ map snd xs)


diffusionMap :: FilePath
             -> FilePath
             -> IO ()
diffusionMap output input = shelly $ run_ "sc_utils"
    ["reduce", "--method", "dm", T.pack input, T.pack output]