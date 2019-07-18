import scrublet as scr
import scipy as sp
import numpy as np

from .Utils import readMatrix

def detectDoublet(args):
    counts_matrix = readMatrix(args.input)
    scrub = scr.Scrublet(counts_matrix, expected_doublet_rate = 0.06)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(
        min_counts=1, 
        min_cells=3, 
        min_gene_variability_pctl=85,
        n_prin_comps=30
        )

    with open(args.output, 'w') as fl:
        fl.write(str(scrub.threshold_))
        fl.write("\n")
        fl.write('\t'.join(map(str, (doublet_scores.tolist()))))
        fl.write("\n")
        fl.write('\t'.join(map(str, scrub.doublet_scores_sim_)))