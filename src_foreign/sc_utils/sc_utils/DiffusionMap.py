import scipy as sp
import numpy as np
import math

from .Utils import InputData, regress

def diffusionMap(args):
    data = InputData(args.input)
    n_dim = 15
    print("Read Data")
    indptr = [0]
    indices = []
    mat = []
    for row in iter(data):
        for (i,x) in row:
            indices.append(i)
            mat.append(1)
        indptr.append(len(indices))
    mat = sp.sparse.csr_matrix((mat, indices, indptr), dtype=int)
    (n,m) = mat.get_shape()

    coverage = mat.sum(axis=1)

    print("Compute similarity matrix")
    jm = mat.dot(mat.T).todense()
    s = coverage.dot(np.ones((1,n)))
    jm = jm / (s + s.T - jm)

    # regression
    K = regression(jm, coverage/m)

    # Gaussian kernel
    #jm = np.exp(- (1/jm) / 0.2)
    #np.fill_diagonal(jm, 0)

    print("Normalization")
    Q_i = np.diag(K.sum(axis=1)**-0.5)
    T_ = np.matmul(np.matmul(Q_i, K), Q_i)

    print("Reduction")
    (evals, evecs) = sp.sparse.linalg.eigs(T_, k=n_dim+1, which='LR')
    ix = evals.argsort()[::-1]
    evals = np.real(evals[ix])
    evecs = np.matmul(Q_i, np.real(evecs[:, ix]))
    dmap = np.matmul(evecs, np.diag(evals))
    np.savetxt(args.output, dmap[:, 1:], delimiter='\t')

def regression(mat, coverage):
    n, m = mat.shape

    X = 1 / coverage.dot(np.ones((1,n)))
    X = 1 / (X + X.T - 1)
    X = X[np.triu_indices(n, k = 1)].T

    y = mat[np.triu_indices(n, k = 1)].T

    print(sp.stats.spearmanr(X[...,0], y[...,0]))
    res = np.zeros((n,m))
    y = y / regress(X,y)
    res[np.triu_indices(n, k = 1)] = y.flatten()
    res = res + res.T
    return res