from sklearn.linear_model import LogisticRegression
from sklearn.metrics import log_loss
import numpy as np
from scipy.stats import chi2

from statsmodels.stats.multitest import multipletests


from .Utils import InputData, readMatrix

def diff(args):
    fg = readMatrix(args.fg)
    bg = readMatrix(args.bg)
    idx, probs = diffTest(fg, bg)
    _,ps,_,_ = multipletests(np.exp(np.array(probs)), 0.1, method="fdr_bh")
    print(list(filter(lambda x: x < 0.01, ps)))

def diffTest(fg, bg):
    (n1,_) = fg.shape
    (n2,_) = bg.shape

    idx1 = np.where(np.sum(fg, axis=0) > 0.1*n1)[1]
    idx2 = np.where(np.sum(bg, axis=0) > 0.1*n2)[1]
    idx = list(set(list(idx1) + list(idx2)))
    idx.sort()

    fg_depth = np.log10(np.sum(fg, axis=1))
    bg_depth = np.log10(np.sum(bg, axis=1))
    z = np.concatenate((fg_depth, bg_depth))
    X = np.array([1] * n1 + [0] * n2)[:, np.newaxis]
    probs = []
    print(len(idx))
    for i in idx:
        if (i % 1000 == 0): print(i)
        a1 = fg[:, i].todense()
        a2 = bg[:, i].todense()
        Y = np.ravel(np.concatenate((a1,a2)))
        probs.append(test(X,Y,z))
    return (idx, probs)

"""
X: sample, 1d vector
Y: label, 1d vector
z: depth
"""
def test(X, Y, z):
    prop_accessible = np.full(X.shape, np.sum(Y) / len(Y))

    X_full = np.concatenate((X, prop_accessible, z), axis=1)
    full = LogisticRegression(penalty="none", random_state=0,
        solver="lbfgs", multi_class='ovr'
        ).fit(X_full, Y)
    y_prob = full.predict_proba(X_full)
    full_prob = log_loss(Y, y_prob, normalize=False)

    X_reduced = np.concatenate((prop_accessible, z), axis=1)
    reduced = LogisticRegression(penalty="none", random_state=0,
        solver="lbfgs", multi_class='ovr'
        ).fit(X_reduced, Y)
    y_prob = reduced.predict_proba(X_reduced)
    reduced_prob = log_loss(Y, y_prob, normalize=False)

    chi = 2 * (reduced_prob - full_prob)
    return chi2.logsf(chi, 1)