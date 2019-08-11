import scipy as sp
import numpy as np
import math
from sklearn.linear_model import LinearRegression

from .Utils import readMatrix, regress

def diffusionMap(args):
    nSample = 2000
    nChunk = 1000

    print("Read Data")
    mat = readMatrix(args.input, binary=True)

    (n,_) = mat.get_shape()
    if (nSample < n):
        idx = np.arange(n)
        np.random.shuffle(idx)
        sample = mat[idx[:nSample], :]
        dm = DiffusionMap(sample, sampling_rate=nSample/n)
        i = nSample
        res = [dm.coordinates]
        while i < n:
            data = mat[idx[i:i+nChunk], :]
            res.append(dm.fit(data))
            i = i + nChunk
        res = np.concatenate(res, axis=0)[:, 1:]
        res = res[np.argsort(idx), :]
    else:
        res = DiffusionMap(mat).coordinates[:, 1:]

    np.savetxt(args.output, res, delimiter='\t')

class DiffusionMap:
    def __init__(self, mat, n_dim=30, sampling_rate=1):
        self.sample = mat
        self.sampling_rate = sampling_rate
        self.dim = mat.get_shape()[1]
        self.coverage = mat.sum(axis=1) / self.dim

        print("Compute similarity matrix")
        jm = jaccard(mat)

        self.normalizer = Normalizer(jm, self.coverage)

        S = self.normalizer.fit(jm, self.coverage, self.coverage)
        np.fill_diagonal(S, 0)

        print("Normalization")
        self.D = np.diag(1/(self.sampling_rate * S.sum(axis=1)))
        L = np.matmul(self.D, S)

        print("Reduction")
        (evals, evecs) = sp.sparse.linalg.eigs(L, n_dim+1, which='LR')
        ix = evals.argsort()[::-1]
        self.evals = np.real(evals[ix])
        self.evecs = np.real(evecs[:, ix])
        self.coordinates = self.evecs
        #self.coordinates = np.matmul(np.matmul(self.Q, self.evecs), np.diag(self.evals**self.t))
        print(self.evals)
        print(self.evecs)
        #self.evecs = np.matmul(Q_i, np.real(evecs[:, ix]))

    def fit(self, data):
        jm = jaccard2(self.sample, data)
        S_ = self.normalizer.fit(jm, self.coverage, data.sum(axis=1) / self.dim).T
        D_ = np.diag(1/(self.sampling_rate * S_.sum(axis=1)))
        L_ = np.matmul(D_, S_)
        evecs = (L_.dot(self.evecs)).dot(np.diag(1/self.evals))
        return evecs

""" Compute pair-wise jaccard index
"""
def jaccard(mat):
    (n,_) = mat.get_shape()
    coverage = mat.sum(axis=1)
    jm = mat.dot(mat.T).todense()
    c = coverage.dot(np.ones((1,n)))
    jm = jm / (c + c.T - jm)
    # Gaussian kernel
    #jm = np.exp(- (1/jm) / 0.2)
    return jm

""" Compute pair-wise jaccard index between two matrix.
Input:
    mat1: n1 x m
    mat2: n2 x m
Output:
    jm: n1 x n2
"""
def jaccard2(mat1, mat2):
    coverage1 = mat1.sum(axis=1)
    coverage2 = mat2.sum(axis=1)

    jm = mat1.dot(mat2.T).todense()
    n1, n2 = jm.shape
    c1 = coverage1.dot(np.ones((1,n2)))
    c2 = coverage2.dot(np.ones((1,n1)))
    jm = jm / (c1 + c2.T - jm)
    return jm

class Normalizer:
    def __init__(self, jm, c):
        n, _ = jm.shape

        X = 1 / c.dot(np.ones((1,n)))
        X = 1 / (X + X.T - 1)
        X = X[np.triu_indices(n, k = 1)].T
        y = jm[np.triu_indices(n, k = 1)].T
        print(sp.stats.spearmanr(X[...,0], y[...,0]))

        self.model = LinearRegression().fit(X, y)

    def fit(self, jm, c1, c2):
        X1 = 1 / c1.dot(np.ones((1, c2.shape[1]))) 
        X2 = 1 / c2.dot(np.ones((1, c1.shape[1])))
        X = 1 / (X1 + X2.T - 1)

        y = self.model.predict(X.flatten().T).reshape(jm.shape)
        return np.array(jm / y)

def scale(mat):
    def scaling(xs):
        s = 0
        for x in xs:
            s = s + x * x
        s = math.sqrt(s)
        return np.array([x / s for x in xs])
    return np.apply_along_axis(scaling, 1, mat)