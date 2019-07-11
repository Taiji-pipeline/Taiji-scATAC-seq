import scipy as sp
import numpy as np
import math

from .Utils import InputData, regressOut

def diffusionMap(args):
    data = InputData(args.input)
    n_dim = 15
    t = 1
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
    (n,_) = mat.get_shape()

    coverage = mat.sum(axis=1)

    print("Compute similarity matrix")
    jm = mat.dot(mat.T).todense()
    s = coverage.dot(np.ones((1,n)))
    jm = jm / (s + s.T - jm)

    print(jm)

    jm = regression(jm, coverage)
    jm = jm.clip(min=0)
    jm = (jm + jm.T) / 2
    print(jm)

    np.fill_diagonal(jm, 0)

    print("Normalization")
    s = jm.sum(axis=1).dot(np.ones((1,n)))
    jm = jm / s

    print("Reduction")
    (evals, evecs) = sp.sparse.linalg.eigs(jm, k=n_dim+1, which='LR')
    ix = evals.argsort()[::-1][1:]
    evals = np.real(evals[ix])
    evecs = np.real(evecs[:, ix])
    dmap = np.matmul(evecs, np.diag(evals**t))
    np.savetxt(args.output, dmap, delimiter='\t')

def regression(mat, coverage):
    s = np.sqrt(np.matmul(coverage, coverage.T))
    X = s.flatten().reshape(s.size, 1)
    y = mat.flatten().reshape(mat.size, 1)
    return regressOut(X,y).reshape(mat.shape)
