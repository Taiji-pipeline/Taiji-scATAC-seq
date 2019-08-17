import scrublet as scr
import scipy as sp
import numpy as np
from sklearn.mixture import BayesianGaussianMixture

import matplotlib.pyplot as plt  

from .Utils import readMatrix

def detectDoublet(args):
    counts_matrix = readMatrix(args.input, binary=False)
    scrub = scr.Scrublet(counts_matrix,
        expected_doublet_rate = 0.06, sim_doublet_ratio=3, n_neighbors=25)
    doublet_scores, _ = scrub.scrub_doublets(
        min_counts=1, 
        min_cells=3, 
        min_gene_variability_pctl=85,
        n_prin_comps=30
        )

    threshold = findThresHold_fit(scrub.doublet_scores_sim_)

    with open(args.output, 'w') as fl:
        fl.write(str(threshold))
        fl.write("\n")
        fl.write('\t'.join(map(str, (doublet_scores.tolist()))))
        fl.write("\n")
        fl.write('\t'.join(map(str, scrub.doublet_scores_sim_)))

def findThresHold_fit(X):
    X = np.array([X]).T
    gmm = BayesianGaussianMixture(n_components=2, max_iter=1000).fit(X)
    i = np.argmax(gmm.means_)
    probs = gmm.predict_proba(X)[:,i]
    vals = X[np.argwhere(probs>0.5)].flatten()
    th = min(vals)
    return th