from sklearn.linear_model import LogisticRegression
from sklearn.metrics import log_loss
import numpy as np
from scipy.stats import chi2
from statsmodels.stats.multitest import fdrcorrection

from .Utils import InputData, readMatrix

def diff(args):
    fg = readMatrix(args.fg, binary=True)
    bg = readMatrix(args.bg, binary=True)
    idx, pval, enrichment = diffTest(fg, bg)
    fdr = fdrcorrection(pval)[1]
    np.savetxt( args.output, np.column_stack((idx, enrichment, pval, fdr)),
        delimiter='\t',
        fmt='%i %10.5f %1.4e %1.4e' )

def process(fg, bg, X, z, idx):
    probs = []
    enrichment = []
    pseudoCount_fg = 1 / fg.shape[0]
    pseudoCount_bg = 1 / bg.shape[0]
    for i in range(len(idx)):
        if (i % 500 == 0): print(i)
        i = idx[i]
        a1 = fg[:, i].todense()
        a2 = bg[:, i].todense()
        Y = np.ravel(np.concatenate((a1,a2)))
        probs.append(test(X,Y,z))
        e = (np.average(a1) + pseudoCount_fg) / (np.average(a2) + pseudoCount_bg)
        enrichment.append(e)
    return (probs, enrichment)

def diffTest(fg, bg):
    (n1,_) = fg.shape
    (n2,_) = bg.shape

    idx1 = np.where(np.sum(fg, axis=0) > 0.05*n1)[1]
    idx2 = np.where(np.sum(bg, axis=0) > 0.05*n2)[1]
    idx = list(set(list(idx1) + list(idx2)))
    idx.sort()

    fg_depth = np.log10(np.sum(fg, axis=1))
    bg_depth = np.log10(np.sum(bg, axis=1))
    z = np.concatenate((fg_depth, bg_depth))
    X = np.array([1] * n1 + [0] * n2)[:, np.newaxis]
    print(len(idx))
    (probs, enrichment) = process(fg, bg, X, z, idx)
    return (idx, probs, enrichment)

"""
X: sample, 1d vector
Y: label, 1d vector
z: depth
"""
def test(X, Y, z):
    prop_accessible = np.full(X.shape, np.sum(Y) / len(Y))

    model = LogisticRegression(penalty="none", random_state=0,
        solver="lbfgs", multi_class='ovr', tol=1e-3, warm_start=False
        )

    X_reduced = np.concatenate((prop_accessible, z), axis=1)
    model.fit(X_reduced, Y)
    y_prob = model.predict_proba(X_reduced)
    reduced_prob = log_loss(Y, y_prob, normalize=False)

    X_full = np.concatenate((X, prop_accessible, z), axis=1)
    model.fit(X_full, Y)
    y_prob = model.predict_proba(X_full)
    full_prob = log_loss(Y, y_prob, normalize=False)

    chi = 2 * (reduced_prob - full_prob)
    return chi2.sf(chi, 1)