{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE DataKinds #-}
module Taiji.Pipeline.SC.ATACSeq.Functions
    ( module Taiji.Pipeline.SC.ATACSeq.Functions.Preprocess
    ) where

import Bio.Utils.BitVector
import Bio.Data.Experiment
import Bio.Data.Experiment.Parser
import           Bio.Pipeline.Utils
import Control.Lens
import Control.Monad.Reader (asks)
import Text.Printf (printf)
import           Bio.HTS
import Conduit
import qualified Data.Text as T
import qualified Data.ByteString.Char8 as B
import Bio.Data.Bed.Utils
import Scientific.Workflow

import Taiji.Pipeline.SC.ATACSeq.Functions.Preprocess
import Taiji.Pipeline.SC.ATACSeq.Types


-- | Remove duplicates for reads originated from a single cell.
rmDup :: BAMHeader -> [BAM] -> [BAM]
rmDup = undefined

rmDuplicates :: SCATACSeqConfig config
             => ATACSeq S (Either (File '[NameSorted] 'Bam)
                         (File '[NameSorted, PairedEnd] 'Bam))
             -> WorkflowConfig config (ATACSeq S (File '[] 'Bam))
rmDuplicates input = do
    dir <- asks _scatacseq_output_dir >>= getPath . (<> (asDir "/Bam"))
    let output = printf "%s/%s_rep%d_filt_dedup.bam" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
        f fl = do
            header <- getBamHeader fl
            runResourceT $ runConduit $ streamBam fl .| getReadGroup .|
                concatMapC (rmDup header) .| sinkBam output header
            return $ location .~ output $ emptyFile
    input & replicates.traverse.files %%~ liftIO . f .
        either (^.location) (^.location)

getReadGroup:: Monad m => ConduitT BAM [BAM] m ()
getReadGroup = go []
  where
    go acc = await >>= \case 
        Nothing -> yield acc
        Just x -> case acc of
            [] -> go [x]
            a:_ -> if getBarcode a == getBarcode x
                then go $ x : acc
                else yield acc >> go [a]
    getBarcode bam = Just $ head $ B.split ':' $ queryName bam


{-
scAtacMkBaseMap :: ATACSeqConfig config
                => ATACSeq S (File '[] 'Bam)
                -> WorkflowConfig config (ATACSeq S (File '[] 'Bam))

mkBaseMap :: FilePath   -- ^ Bam file
          -> [(B.ByteString, Int)]
          -> Int  -- ^ resolution
          -> IO (M.HashMap B.ByteString (U.Vector Bool))
mkBaseMap input chrSize res = runResourceT (
    runConduit $ streamBam input .| bamToBedC .| baseMap chrSize) >>=
    (\BaseMap x -> fmap (downSample res) x) 

-- | Down-sample basemap to given resolution by max-pooling.
downSample :: Int    -- ^ Resolution, must be the multiple of 8.
           -> BitVector
           -> U.Vector Bool
downSample res (BitVector _ v) = U.generate n $ \i ->
    U.any (/=0) $ U.slice (i * step) step v
  where
    step = res `div` 8
    n = if U.length v `mod` step == 0
        then U.length v `div` step
        else U.length v `div` step + 1
        -}