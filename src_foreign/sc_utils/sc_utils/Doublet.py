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

    threshold = findThresHold_fit(scrub.doublet_scores_sim_, output=args.plot)

    with open(args.output, 'w') as fl:
        fl.write(str(threshold))
        fl.write("\n")
        fl.write('\t'.join(map(str, (doublet_scores.tolist()))))
        fl.write("\n")
        fl.write('\t'.join(map(str, scrub.doublet_scores_sim_)))

def findThresHold_fit(X, output=None):
    X = np.array([X]).T
    gmm = BayesianGaussianMixture(n_components=2).fit(X)
    
    m1 = gmm.means_[0]
    s1 = gmm.covariances_[0,0]
    m2 = gmm.means_[1]
    s2 = gmm.covariances_[1,0]

    c1 = sp.stats.norm.ppf(0.9,m1,s1)[0]
    c2 = sp.stats.norm.ppf(0.01,m2,s2)[0]
    th = max(c1,c2)

    if (output):
        plt.hist(X, bins=150, density=True)
        drawGaussian(m1,s1)
        drawGaussian(m2,s2)
        plt.axvline(x=th)
        plt.savefig(output)

    return th

def drawGaussian(m, s):
    x = np.linspace(0, 1, 150)
    y = sp.stats.norm.pdf(x,m,s)
    plt.plot(x,y)