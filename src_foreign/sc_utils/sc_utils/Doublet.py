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
    probs = gmm.predict_proba(X)[:,1]
    th = min(X[np.argwhere(probs>0.9)])[0][0]
    
    '''
    m1 = gmm.means_[0]
    s1 = gmm.covariances_[0,0]
    m2 = gmm.means_[1]
    s2 = gmm.covariances_[1,0]

    c1 = sp.stats.norm.ppf(0.9,m1,s1)[0]
    c2 = sp.stats.norm.ppf(0.01,m2,s2)[0]
    th = max(c1,c2)
    print(th)
    '''

    return th