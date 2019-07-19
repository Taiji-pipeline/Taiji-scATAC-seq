import scrublet as scr
import scipy as sp
import numpy as np
from sklearn.neighbors.kde import KernelDensity

from .Utils import readMatrix

def detectDoublet(args):
    counts_matrix = readMatrix(args.input)
    scrub = scr.Scrublet(counts_matrix,
        expected_doublet_rate = 0.06, sim_doublet_ratio=3, n_neighbors=25)
    doublet_scores, _ = scrub.scrub_doublets(
        min_counts=1, 
        min_cells=3, 
        min_gene_variability_pctl=85,
        n_prin_comps=30
        )
    threshold = findThresHold(scrub.doublet_scores_sim_)

    with open(args.output, 'w') as fl:
        fl.write(str(threshold))
        fl.write("\n")
        fl.write('\t'.join(map(str, (doublet_scores.tolist()))))
        fl.write("\n")
        fl.write('\t'.join(map(str, scrub.doublet_scores_sim_)))

def findThresHold(X):
    X = np.array([X]).T
    kde = KernelDensity(kernel='gaussian', bandwidth=0.005).fit(X)
    sample_X = np.linspace(0, 1, 500)[:, np.newaxis]
    density = np.exp(kde.score_samples(sample_X))
    gradient = np.gradient(density)
    return sample_X[np.where(gradient < 0.0001)[0][0],0]