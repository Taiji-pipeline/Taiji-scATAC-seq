import scipy as sp
import numpy as np
import math

from .Utils import InputData, regress

def diffusionMap(args):
    data = InputData(args.input)
    n_dim = 15
    t = 5
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

    jm = regression(jm, coverage/m)
    jm = (jm + jm.T) / 2
    jm = jm.clip(min=0)

    #np.fill_diagonal(jm, 0)

    print("Normalization")
    s = jm.sum(axis=1).dot(np.ones((1,n)))
    jm = jm / s

    print("Reduction")
    (evals, evecs) = sp.sparse.linalg.eigs(jm, k=n_dim+1, which='LR')
    ix = evals.argsort()[::-1]
    evals = np.real(evals[ix])
    evecs = np.real(evecs[:, ix])
    dmap = np.matmul(evecs, np.diag(evals**t))
    print(dmap[:,0])
    np.savetxt(args.output, dmap[:, 1:], delimiter='\t')

def regression(mat, coverage):
    n, m = mat.shape

    X = 1 / coverage.dot(np.ones((1,n)))
    X = 1 / (X + X.T - 1)
    X = X.reshape(m*n, 1)

    y = mat.flatten().reshape(m*n, 1)

    print(sp.stats.spearmanr(X, y))
    y_ = y / regress(X,y)
    print(sp.stats.spearmanr(X, y_.reshape(m*n, 1)))
    return y_.reshape(n,m)

