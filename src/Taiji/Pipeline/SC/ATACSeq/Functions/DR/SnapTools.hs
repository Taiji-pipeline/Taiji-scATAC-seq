{-# LANGUAGE DataKinds #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE QuasiQuotes       #-}
{-# LANGUAGE TemplateHaskell   #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE ViewPatterns #-}
module Taiji.Pipeline.SC.ATACSeq.Functions.DR.SnapTools
   ( snapPre
   , performSnap
   , mkSnapMat
   ) where

import qualified Data.ByteString.Char8 as B
import Bio.Seq.IO
import Bio.Utils.Misc (readInt)
import Data.ByteString.Lex.Integral (packDecimal)
import Bio.Data.Bed
import Shelly hiding (FilePath)
import qualified Data.Text as T
import System.IO.Temp (withTempFile)
import System.IO
import Data.Conduit.Internal (zipSources)
import qualified Data.Vector as V

import qualified Language.R                        as R
import           Language.R.QQ
import Language.R.HExp

import Taiji.Prelude
import Taiji.Pipeline.SC.ATACSeq.Types
import Taiji.Pipeline.SC.ATACSeq.Functions.Utils

snapPre :: SCATACSeqConfig config
        => SCATACSeq S (File '[NameSorted, Gzip] 'Bed)
        -> ReaderT config IO (SCATACSeq S (File '[] 'Other))
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
                , "--bin-size-list", "5000" ]
            return $ location .~ output $ emptyFile )

performSnap :: SCATACSeqConfig config 
            => SCATACSeq S (File '[] 'Other)
            -> ReaderT config IO (SCATACSeq S (File '[] 'Tsv, File '[Gzip] 'Tsv))
performSnap input = input & replicates.traverse.files %%~ ( \fl -> do
    dir <- asks ((<> "/Snap/") . _scatacseq_output_dir) >>= getPath
    let output1 = printf "%s/%s_rep%d_snap_rownames.tsv" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
        output2 = printf "%s/%s_rep%d_snap.tsv.gz" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    liftIO $ snap output1 output2 (fl^.location) 
    return ( location .~ output1 $ emptyFile
           , location .~ output2 $ emptyFile )
    )

snap :: FilePath   -- ^ Row names
     -> FilePath   -- ^ Matrix
     -> FilePath   -- ^ Input
     -> IO ()
snap rownames mat input = R.runRegion $ do
    _ <- [r| library("SnapATAC")
        x.sp <- createSnap(file=input_hs, sample="SNAP")
        x.sp <- addBmatToSnap(x.sp, bin.size=5000, num.cores=1)
        x.sp <- makeBinary(x.sp, mat="bmat")
        x.sp <- runJDA( obj=x.sp, input.mat="bmat",
            bin.cov.zscore.lower=-2,
            bin.cov.zscore.upper=2,
            pc.num=30,
            norm.method="normOVE",
            max.var=5000,
            do.par=TRUE,
            ncell.chunk=1000,
            num.cores=10,
            seed.use=10,
            tmp.folder=tempdir()
        )
        write.table(cbind(x.sp@barcode, rowSums(x.sp@bmat)),
            file=rownames_hs, sep="\t", row.names=F, col.names=F, quote=F)
        gz1 <- gzfile(mat_hs, "w")
        write.table(x.sp@smat@dmat, gz1, sep="\t", row.names=F, col.names=F, quote=F)
        close(gz1)
    |]
    return ()


mkSnapMat :: SCATACSeqConfig config 
          => SCATACSeq S (File t 'Other)
          -> ReaderT config IO (SCATACSeq S (File '[] 'Tsv, File '[Gzip] 'Tsv))
mkSnapMat input = input & replicates.traverse.files %%~ ( \x -> do
    dir <- asks ((<> "/Snap/") . _scatacseq_output_dir) >>= getPath
    let output1 = printf "%s/%s_rep%d_snap_rownames_new.tsv" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
        output2 = printf "%s/%s_rep%d_snap_new.tsv.gz" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    liftIO $ snap' output1 output2 (x^.location)
    return ( location .~ output1 $ emptyFile
           , location .~ output2 $ emptyFile )
    )


snap' :: FilePath -- ^ Row names
      -> FilePath -- ^ Mat output
      -> FilePath   -- ^ Sparse Matrix
      -> IO ()
snap' rownames matOutput matFl = withTempDir (Just "./") $ \tmpdir -> do
    --regions <- runResourceT $ runConduit $ streamBedGzip nameFl .| mapC f .| sinkList
    let ridx =  tmpdir <> "/ridx.txt"
        cidx =  tmpdir <> "/cidx.txt"
        vals = tmpdir <> "/vals.txt"
        rname = tmpdir <> "/rname.txt"
    mat <- mkSpMatrix readInt matFl
    parseSpMat rname ridx cidx vals mat
    putStrLn "Run Snap..."
    R.runRegion $ do
        _ <- [r| library("SnapATAC")
                 library("Rcpp")
                 library("inline")
                 i <- scan(ridx_hs, integer())
                 j <- scan(cidx_hs, integer())
                 vals <- scan(vals_hs, integer())
                 M <- sparseMatrix(i,j,x=vals)
                 rm(i)
                 rm(j)
                 rm(vals)
                 rn <- scan(rname_hs, character())
                 rownames(M) <- rn
                 saveRDS(M, file="mat.rds")
                 M <- t(M)

                 sparseProdCpp <- '
                     using Eigen::MappedSparseMatrix;
                     using Eigen::SparseMatrix;
                     const MappedSparseMatrix<double>  A(as<MappedSparseMatrix<double> >(AA));
                     const MappedSparseMatrix<double>  B(as<MappedSparseMatrix<double> >(BB));
                     return wrap(A * B.adjoint());
                     '
                 incl <- '
                     using   Eigen::LLT;
                     using   Eigen::Lower;
                     using   Eigen::Map;
                     using   Eigen::MatrixXd;
                     using   Eigen::MatrixXi;
                     using   Eigen::Upper;
                     using   Eigen::VectorXd;
                     typedef Map<MatrixXd>  MapMatd;
                     typedef Map<MatrixXi>  MapMati;
                     typedef Map<VectorXd>  MapVecd;
                     inline  MatrixXd AtA(const MatrixXd& A) {
                         int    n(A.cols());
                         return   MatrixXd(n,n).setZero().selfadjointView<Lower>()
                                 .rankUpdate(A.adjoint());
                     }
                     inline  MatrixXd AAt(const MatrixXd& A) {
                         int    n(A.cols());
                         return   MatrixXd(n,n).setZero().selfadjointView<Lower>()
                                 .rankUpdate(A);
                     }
                     '
                #sparseProd <- cxxfunction(signature(AA = "dgCMatrix", BB = "dgCMatrix"), sparseProdCpp, "RcppEigen", incl)

                runJaccard2 <- function(obj1, obj2 ) {
                    calJaccard <- function(X_i, X_j){
                        #A = sparseProd(X_i, X_j);
                        A = Matrix::tcrossprod(X_i, X_j);
                        bi = Matrix::rowSums(X_i);
                        bj = Matrix::rowSums(X_j);
                        jmat = as.matrix(A / (replicate(ncol(A), bi) + t(replicate(nrow(A), bj)) - A));
                        rm(A);
                        rm(bi);
                        rm(bj);
                        gc();	
                        return(jmat);				
                    }
	
                    jmat = calJaccard(obj1@bmat, obj2@bmat);

                    # remove large objects
                    obj1@jmat@jmat = jmat;
                    obj1@jmat@p1 = Matrix::rowMeans(obj1@bmat);
                    obj1@jmat@p2 = Matrix::rowMeans(obj2@bmat);
                    obj1@jmat@norm = FALSE;
                    obj1@jmat@method = character();
                    rm(jmat);
                    gc();
                    return(obj1);
                }

                eig_decomp <- function(M, n_eigs, sym = isSymmetric(M)) {
                    n <- nrow(M)
                    f <- function(x, A = NULL) as.matrix(A %*% x)
                    wh <- if (sym) 'LA' else 'LM'
                    #constraints: n >= ncv > nev
                    ar <- igraph::arpack(f, extra = M, sym = sym, options = list(
                        which = wh, n = n, ncv = min(n, 4*n_eigs), nev = n_eigs + 1))
                    if (!sym) {
                        ar$vectors <- Re(ar$vectors)
                        ar$values  <- Re(ar$values)
                    }
                    if (length(dim(ar$vectors)) == 0L)
                        ar$vectors <- matrix(ar$vectors, ncol = 1L)
                    ar
                }

                .normOVE <- function(p1, p2){
                    pp = tcrossprod(p1, p2);
                    ss = matrix(rep(p1,each=length(p2)), ncol=length(p2), byrow=TRUE) +  matrix(rep(p2, each=length(p1)), ncol=length(p2), byrow=FALSE)
                    ee = pp/(ss - pp)
                    return(ee)	
                }

                num_landmark = 20000
                seed_use = 3944

                # create snap object
                x.sp = newSnap();
                x.sp@barcode = colnames(M);
                x.sp@bmat = as(t(M), "dgCMatrix");
                print(nrow(x.sp))
                if(max(x.sp@bmat) > 1){
                    x.sp = makeBinary(x.sp, mat="bmat");		
                }

                bin.cov = log(Matrix::colSums(x.sp@bmat)+1,10);
                x.sp = x.sp[,which(bin.cov > 0), mat="bmat"];
                bin.cov = bin.cov[which(bin.cov > 0)]
                bin.covs.zscore = (bin.cov - mean(bin.cov)) / sd(bin.cov);
                cutoff = sort(bin.covs.zscore)[length(bin.covs.zscore) * 0.95];
                x.sp = x.sp[,which(bin.covs.zscore <= cutoff), mat="bmat"];

                # 2. Randomly Sample 10k Cells as Landmark cells
                row.covs = log((Matrix::rowSums(x.sp@bmat))+1,10);		
                row.covs.dens <- density(x = row.covs, bw = 'nrd', adjust = 1)
                sampling_prob <- 1 / (approx(x = row.covs.dens$x, y = row.covs.dens$y, xout = row.covs)$y + .Machine$double.eps)
                set.seed(as.numeric(seed_use))
                idx.landmark.ds <- sort(sample(x = seq(row.covs), size = num_landmark, prob = sampling_prob));
                x.landmark.sp = x.sp[idx.landmark.ds,]
                x.landmark.sp = runJaccard2(x.landmark.sp, x.landmark.sp);

                # 3 train a regression model
                row.covs = log((Matrix::rowSums(x.landmark.sp@bmat))+1,10);		
                row.covs.dens <- density(x = row.covs, bw = 'nrd', adjust = 1)
                sampling_prob <- 1 / (approx(x = row.covs.dens$x, y = row.covs.dens$y, xout = row.covs)$y + .Machine$double.eps)
                set.seed(as.numeric(seed_use))
                idx.ds <- sort(sample(x = seq(row.covs), size = min(2000, num_landmark), prob = sampling_prob));

                jmat.tr = x.landmark.sp@jmat@jmat[idx.ds,idx.ds];
                b1.tr = x.landmark.sp@jmat@p1[idx.ds];
                b2.tr = x.landmark.sp@jmat@p2[idx.ds];
                emat.tr = .normOVE(b1.tr, b2.tr);
                data = data.frame(x=emat.tr[upper.tri(emat.tr)], y=jmat.tr[upper.tri(jmat.tr)])	
                model <- lm(y ~ x + I(x^2), data);

                # 6 Perform diffuion map against landmark cells
                b1.te = x.landmark.sp@jmat@p1;
                b2.te = x.landmark.sp@jmat@p2;
                jmat.te = x.landmark.sp@jmat@jmat;
                emat.te = .normOVE(b1.te, b2.te);

                preds = predict(model, data.frame(x=emat.te[upper.tri(emat.te)]), se.fit = TRUE);
                norm = as.numeric(jmat.te[upper.tri(jmat.te)]/preds$fit);
                nmat = matrix(0, nrow(jmat.te), ncol(jmat.te));
                nmat[upper.tri(nmat)]  = norm;
                nmat[lower.tri(nmat)] = t(nmat)[lower.tri(nmat)];

                diag(nmat) = 0
                norm_p1 = nmat
                d_norm1 <- Matrix::rowSums(norm_p1);
                d_rot1 <- Diagonal(x = d_norm1 ^ -.5);
                transitions <- as.matrix(d_rot1 %*% norm_p1 %*% d_rot1);
                diag(transitions) = 0;
                eig_transitions = eig_decomp(transitions, 50);
                x.landmark.sp@smat@dmat = eig_transitions$vectors;
                x.landmark.sp@smat@sdev = eig_transitions$value;


                x.os.sp = x.sp[-idx.landmark.ds,];
                chunk.size = 10000
                idx.ls = split(seq(nrow(x.os.sp)), ceiling(seq_along(seq(nrow(x.os.sp)))/chunk.size))
                eigen.vector.mat = x.landmark.sp@smat@dmat;

                eigen.vector.mat.list = list()
                for(i in seq(idx.ls)){
                    idx = idx.ls[[i]];
                    print(paste("chunk", i, sep=" "))
                    x.sub.sp = x.os.sp[idx,];
                    x.sub.sp = runJaccard2(x.sub.sp, x.landmark.sp);
                    
                    b1.te = x.sub.sp@jmat@p1;
                    b2.te = x.sub.sp@jmat@p2;
                    jmat.te = x.sub.sp@jmat@jmat;
                    emat.te = .normOVE(b1.te, b2.te);

                    preds = predict(model, data.frame(x=c(emat.te)), se.fit = TRUE);
                    norm = as.numeric(c(jmat.te)/preds$fit);
                    nmat = matrix(norm, nrow(jmat.te), ncol(jmat.te))

                    trans_p = nmat;
                    norm_p = trans_p
                    d_new <- Matrix::rowSums(trans_p, na.rm = TRUE)

                    d_norm_new <- Matrix::rowSums(norm_p)
                    d_rot_new <- Diagonal(x = d_norm_new ^ -.5)
                    M_new <- d_rot_new %*% norm_p %*% d_rot1
                    eigen.vector = t(t(M_new %*% x.landmark.sp@smat@dmat) / x.landmark.sp@smat@sdev)
                    #eigen.vector.mat = rbind(eigen.vector.mat, eigen.vector);
                    eigen.vector.mat.list[[i]] = eigen.vector
                }


                x.sp@smat@dmat = matrix(0, nr=nrow(x.sp), nc=50)
                x.sp@smat@dmat[idx.landmark.ds,] = x.landmark.sp@smat@dmat[,2:51]
                x.sp@smat@dmat[-idx.landmark.ds,] = as.matrix(do.call(rbind, eigen.vector.mat.list)[,2:51])
                x.sp@smat@sdev = x.landmark.sp@smat@sdev[2:51]

                write.table(cbind(x.sp@barcode, rowSums(x.sp@bmat)),
                    file=rownames_hs, sep="\t", row.names=F, col.names=F, quote=F)
                gz1 <- gzfile(matOutput_hs, "w")
                write.table(x.sp@smat@dmat, gz1, sep="\t", row.names=F, col.names=F, quote=F)
                close(gz1)
        |]
        return ()

parseSpMat :: FilePath  -- ^ row names
           -> FilePath  -- ^ row index
           -> FilePath  -- ^ col index
           -> FilePath  -- ^ value
           -> SpMatrix Int
           -> IO ()
parseSpMat rname ridx cidx vals mat = do
    _ <- runResourceT $ runConduit $
        zipSources (iterateC succ 0) (streamRows mat) .| mapC f .|
        getZipSink ((,,,) <$> ZipSink sink1 <*> ZipSink sink2 <*> ZipSink sink3 <*> ZipSink sink4)
    return ()
  where
    sink1 = mapC fst .| intersperseC " " .| sinkFile rname
    sink2 = concatMapC (map (fromJust . packDecimal) . (^._1) . snd) .| intersperseC " " .| sinkFile ridx
    sink3 = concatMapC (map (fromJust . packDecimal) . (^._2) . snd) .| intersperseC " " .| sinkFile cidx
    sink4 = concatMapC (map (fromJust . packDecimal) . (^._3) . snd) .| intersperseC " " .| sinkFile vals
    f (idx, (bc, xs)) = (bc
        , unzip3 $ map (\(i, (j,_)) -> (i+1,j+1,1)) $ zip (repeat idx) xs)