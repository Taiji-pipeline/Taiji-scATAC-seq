{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE QuasiQuotes       #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE TemplateHaskell   #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.Cicero
  (cicero) where

import Bio.Data.Bed
import Bio.Seq.IO
import Bio.Utils.Misc (readInt)
import qualified Data.ByteString.Char8 as B
import Data.ByteString.Lex.Integral (packDecimal)
import qualified Data.Vector as V

import qualified Language.R                        as R
import           Language.R.QQ
import Language.R.HExp

import Taiji.Prelude
import Taiji.Pipeline.SC.ATACSeq.Functions.Utils
import Taiji.Pipeline.SC.ATACSeq.Types

cicero :: SCATACSeqConfig config
       => ( Maybe (File '[Gzip] 'NarrowPeak)
          , [SCATACSeq S (File '[Gzip] 'Other)] )
       -> ReaderT config IO (Maybe FilePath)
cicero (Just peakFl, [matFl]) = do
    dir <- asks ((<> "/Cicero") . _scatacseq_output_dir) >>= getPath
    let output = dir <> "/cicero.tsv"
    genome <- asks (fromJust . _scatacseq_genome_index)
    liftIO $ do
        chrSizes <- withGenome genome $ return . getChrSizes
        peaks <- runResourceT $ runConduit $ streamBedGzip (peakFl^.location) .| sinkList
        mat <- mkSpMatrix readInt $ matFl^.replicates._2.files.location
        withTempDir Nothing $ \tmpdir -> do
            --let count = tmpdir <> "/count"
            let count = "count.txt"
                chr = tmpdir <> "/chr"
            B.writeFile chr $ B.unlines $
                map (\(a,b) -> a <> "\t" <> B.pack (show b)) chrSizes
            mkAtacCDS count peaks mat
            getConnection output chr count
    return $ Just output
cicero _ = return Nothing

-- | Create a sparse matrix that can be used as the input for "make_atac_cds"
-- in Cicero.
-- Specifically, the file is a tab-delimited text file with three columns.
-- The first column is the peak coordinates in the form
-- "chr10_100013372_100013596", the second column is the cell name, and the
-- third column is an integer that represents the number of reads from that
-- cell overlapping that peak. The file should not have a header line.
-- For example:
-- chr10_100002625_100002940	cell1	1
-- chr10_100006458_100007593	cell2	2
-- chr10_100006458_100007593	cell3	1
mkAtacCDS :: FilePath
          -> [BED3]         -- ^ Peak list
          -> SpMatrix Int   -- ^ Cell by Peak matrix
          -> IO ()
mkAtacCDS output peaks mat = runResourceT $ runConduit $
    streamRows mat .| concatMapC procRow .| unlinesAsciiC .| sinkFile output
  where
    procRow (cell, xs) = map f xs
      where
        f (i, c) = B.intercalate "\t" [colnames V.! i, cell, fromJust $ packDecimal c]
    colnames = V.fromList $ map f peaks
      where
        f x = B.intercalate "_"
            [x^.chrom, fromJust $ packDecimal $ x^.chromStart, fromJust $ packDecimal $ x^.chromEnd]

-- |
--        ddrtree_coord1	ddrtree_coord2
-- cell1    -0.7084047      -0.7232994
-- cell2    -4.4767964       0.8237284
-- cell3     1.4870098      -0.4723493
getConnection :: FilePath  -- ^ Output
              -> FilePath  -- ^ Chromosome size file
              -> FilePath  -- ^ Count matrix
              -> IO ()
getConnection output chr cds = R.runRegion $ do
    _ <- [r| library(cicero)
        input_cds <- make_atac_cds(cds_hs, binarize = TRUE)
        input_cds <- reduceDimension(input_cds, max_components = 2, num_dim=6,
                      reduction_method = 'tSNE', norm_method = "none")
        tsne_coords <- t(reducedDimA(input_cds))
        row.names(tsne_coords) <- row.names(pData(input_cds))
        cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = tsne_coords)
        conns <- run_cicero(cicero_cds, chr_hs)
        write.table(conns, file=output_hs)
    |]
    return ()