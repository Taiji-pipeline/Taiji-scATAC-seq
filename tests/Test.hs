{-# LANGUAGE DataKinds #-}
{-# LANGUAGE OverloadedLists #-}
{-# LANGUAGE OverloadedStrings #-}
module Test
    ( tests
    ) where

import Conduit
import qualified Data.ByteString.Char8 as B
import           Test.Tasty
import           Data.Conduit.Zlib           (ungzip, multiple)
import           Test.Tasty.HUnit
import Bio.Data.Bed
import Bio.Utils.Misc
import Taiji.Prelude
import Data.Function (on)

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

tests :: TestTree
tests = testGroup "Test"
    [ tsseTest
    ]

tsseTest :: TestTree
tsseTest = testCase "TSSe" $ do
    tss <- readPromoter' "data/gencode.gtf.gz"
    fragments <- fmap (groupBy ((==) `on` (^.name))) $
        runResourceT $ runConduit $ streamBedGzip input .| sinkList
    expected @=? map (tssEnrichment tss) fragments
  where
    expected = [11.857707509881424, 2.727272727272727, 6.583072100313478
        , 2.727272727272727, 0.0, 0.0, 1.8181818181818181, 6.1633281972265
        , 0.9090909090909091, 6.220095693779905, 5.965909090909091
        , 7.204116638078901, 9.312638580931262]
    input = "data/fragments.bed.gz"