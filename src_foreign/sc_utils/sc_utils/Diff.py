from sklearn.linear_model import LogisticRegression
from sklearn.metrics import log_loss
import numpy as np
from scipy.stats import chi2
from statsmodels.stats.multitest import fdrcorrection

from .Utils import InputData, readMatrix

def diff(args):
    fg = readMatrix(args.fg, binary=False)
    bg = readMatrix(args.bg, binary=False)
    if (args.index):
        with open(args.index, 'r') as fl:
            idx_set = set([int(l.strip()) for l in fl])
        idx, pval, enrichment = diffTest(fg, bg, idx_set)
    else:
        idx, pval, enrichment = diffTest(fg, bg)
    fdr = fdrcorrection(pval)[1]
    np.savetxt( args.output, np.column_stack((idx, enrichment, pval, fdr)),
        delimiter='\t',
        fmt='%i %10.5f %1.4e %1.4e' )

def diffTest(fg, bg, idx_set=None):
    (n1,m) = fg.shape
    (n2,_) = bg.shape

    fg_depth = np.sum(fg, axis=1)
    fg_scale = np.mean(fg_depth) / 10000
    bg_depth = np.sum(bg, axis=1)
    bg_scale = np.mean(bg_depth) / 10000

    fg_frac = np.ravel(np.sum(fg, axis=0)) / (n1 * fg_scale)
    bg_frac = (np.ravel(np.sum(bg, axis=0)) + 1) / (n2 * bg_scale)

    z = np.concatenate((fg_depth, bg_depth))
    X = np.array([1] * n1 + [0] * n2)[:, np.newaxis]

    def test(i):
        enrichment = fg_frac[i] / bg_frac[i] 
        if enrichment >= 2 and fg_frac[i] > 0:
            a1 = np.clip(fg[:, i].todense(), 0, 1)
            a2 = np.clip(bg[:, i].todense(), 0, 1)
            Y = np.ravel(np.concatenate((a1,a2)))
            prob = getPvalue(X,Y,z)
            return (i, prob, enrichment)
        else:
            return None 
    
    res = []
    if(idx_set):
        idx = idx_set
    else:
        idx = range(m)
    for i in idx:
        r = test(i)
        if r: res.append(r)
    return list(zip(*res))

"""
X: sample, 1d vector
Y: label, 1d vector
z: depth
"""
def getPvalue(X, Y, z):
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