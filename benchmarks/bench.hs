{-# LANGUAGE OverloadedStrings #-}

import Criterion.Main
import Conduit
import qualified Data.ByteString.Char8 as B
import           Data.Conduit.Zlib           (ungzip, multiple)
import Bio.Data.Bed
import Bio.Utils.Misc
import Data.Conduit.List (groupBy, chunksOf)
import Control.Parallel.Strategies (parMap, rseq)

import Taiji.Prelude hiding (groupBy)
import Taiji.Pipeline.SC.ATACSeq.Functions.QC (tssEnrichment, readPromoter)

readPromoter' :: FilePath -> IO (BEDTree Bool)
readPromoter' fl = fmap (bedToTree const) $ runResourceT $ runConduit $
    sourceFile fl .| multiple ungzip .| linesUnboundedAsciiC .| mapC (f . B.split '\t') .| sinkList
  where
    f xs = let chr = xs !! 0
               left = readInt (xs !! 3) - 1
               right = readInt (xs !! 4) - 1
            in case xs !! 6 of
                "-" -> (BED3 chr (right - 2000) $ right + 2001, False)
                _ -> (BED3 chr (left - 2000) $ left + 2001, True)

tsseBatch :: BEDTree Bool -> [BED] -> [Double]
tsseBatch tss fragments = runIdentity $ runConduit $ yieldMany fragments .|
    groupBy ((==) `on` (^.name)) .| mapC (tssEnrichment tss) .| sinkList

tsseBatch' :: BEDTree Bool -> [BED] -> [Double]
tsseBatch' tss fragments = runIdentity $ runConduit $ yieldMany fragments .|
    groupBy ((==) `on` (^.name)) .| chunksOf 20 .|
    concatMapC (parMap rseq (tssEnrichment tss)) .| sinkList

setupEnv = do
    tss <- readPromoter' "data/gencode.gtf.gz"
    fragments <- runResourceT $ runConduit $ streamBedGzip input .| sinkList
    return (tss, fragments, concat $ replicate 100 fragments)
  where
    input = "data/fragments.bed.gz"

main :: IO ()
main = defaultMain [
    env setupEnv $ \ ~(tss, fragments, fragmentsMany) -> bgroup "main"
        [ --bench "tss enrichment" $ whnf (tssEnrichment tss) fragments
         --bench "tss enrichment (batch)" $ nf (tsseBatch tss) fragmentsMany ]
         bench "tss enrichment (batch) - parallel" $ nf (tsseBatch' tss) fragmentsMany ]
    ]