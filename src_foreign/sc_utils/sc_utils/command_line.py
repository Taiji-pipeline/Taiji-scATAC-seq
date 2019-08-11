import argparse
import sc_utils as sc

from .Doublet import detectDoublet
from .Diff import diff

################################################################################
## ARGUMENT PARSER
################################################################################

parser = argparse.ArgumentParser(description="Fast online algorithms")
subparsers = parser.add_subparsers(title="sub commands")

# create the parser for the "run" command
parser_run = subparsers.add_parser('reduce', help='dimension reduction')
parser_run.add_argument('input', type=str, help='gzipped input file')
parser_run.add_argument('output', type=str, help='output matrix in .npy format')
parser_run.add_argument('--method', help='algorithm: svd, lda, dm')
parser_run.add_argument('--sample-size', default=50000, type=int, help='sampling size')
parser_run.set_defaults(func=sc.reduceDimension)

# create the parser for the "clust" command
parser_clust = subparsers.add_parser('clust', help='perform clustering')
parser_clust.add_argument('input', type=str, help='input matrix in .npy format')
parser_clust.add_argument('output', type=str, help='output file')
parser_clust.add_argument('--coverage', help='coverage file')
parser_clust.add_argument('--embed', help='embedding file')
parser_clust.add_argument('--embed-method', help='embedding method')
parser_clust.add_argument('--discard', action='store_true', help='remove first dimension')
parser_clust.add_argument('--scale', action='store_true', help='scale to unit ball')
parser_clust.add_argument('--dim', type=int, help='dimension')
parser_clust.add_argument('-k', default=20, type=int, help='neighbors')
parser_clust.add_argument('--res', type=float, help='resolution')
parser_clust.set_defaults(func=sc.clustering)

# create the parser for the "doublet" command
parser_doublet = subparsers.add_parser('doublet', help='doublet detection')
parser_doublet.add_argument('input', type=str, help='input matrix')
parser_doublet.add_argument('output', type=str, help='output')
parser_doublet.add_argument('--plot', type=str, help='plot')
parser_doublet.set_defaults(func=detectDoublet)

# create the parser for the "diff" command
parser_diff = subparsers.add_parser('diff', help='diff')
parser_diff.add_argument('--fg', type=str, help='foreground matrix')
parser_diff.add_argument('--bg', type=str, help='background matrix')
parser_diff.add_argument('--index', type=str, help='selected index')
parser_diff.add_argument('output', type=str, help='output')
parser_diff.set_defaults(func=diff)

def main():
    args = parser.parse_args()
    args.func(args)