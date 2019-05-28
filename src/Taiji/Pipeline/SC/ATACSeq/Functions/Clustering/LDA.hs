module Taiji.Pipeline.SC.ATACSeq.Functions.Clustering.LDA where where

import qualified Language.R                        as R
import           Language.R.QQ
import qualified Data.Vector.SEXP as V
import Language.R.HExp

cisTopic :: FilePath -> IO [CellCluster]
cisTopic input = do
    members <- R.runRegion $ do
        _ <- [r| library("cisTopic")
                 counts_mel <-
                 cisTopicObject <- createcisTopicObject(counts_mel)
                 rm(counts_mel)
                 cisTopicObject <- runModels(cisTopicObject,
                    topic=50, seed=987, nCores=4,
                    burnin = 250, iterations = 500, addModels=FALSE)

        |]