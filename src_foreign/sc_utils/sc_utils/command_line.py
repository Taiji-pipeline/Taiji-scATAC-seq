import argparse
import sc_utils as sc

################################################################################
## ARGUMENT PARSER
################################################################################

parser = argparse.ArgumentParser(description="Fast online algorithms")
subparsers = parser.add_subparsers(title="sub commands")

# create the parser for the "run" command
parser_run = subparsers.add_parser('reduce', help='dimension reduction')
parser_run.add_argument('input', type=str, help='gzipped input file')
parser_run.add_argument('output', type=str, help='output matrix in .npy format')
parser_run.add_argument('--method', help='algorithm: svd, lda')
parser_run.set_defaults(func=sc.reduceDimension)

# create the parser for the "clust" command
parser_clust = subparsers.add_parser('clust', help='perform clustering')
parser_clust.add_argument('input', type=str, help='input matrix in .npy format')
parser_clust.add_argument('output', type=str, help='output file')
parser_clust.add_argument('--coverage', help='coverage file')
parser_clust.add_argument('--embed', help='embedding file')
parser_clust.add_argument('--embed-method', help='embedding method')
parser_clust.add_argument('--discard', action='store_true', help='remove first dimension')
parser_clust.set_defaults(func=sc.clustering)

def main():
    args = parser.parse_args()
    args.func(args)