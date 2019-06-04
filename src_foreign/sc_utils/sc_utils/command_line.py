import argparse
import sc_utils as sc

################################################################################
## ARGUMENT PARSER
################################################################################

parser = argparse.ArgumentParser(description="Fast online SVD")
subparsers = parser.add_subparsers(title="sub commands")

# create the parser for the "run" command
parser_run = subparsers.add_parser('svd', help='run SVD')
parser_run.add_argument('input', type=str, help='gzipped input file')
parser_run.add_argument('output', type=str, help='output matrix in .npy format')
parser_run.set_defaults(func=sc.lsiTransform)

# create the parser for the "clust" command
parser_clust = subparsers.add_parser('clust', help='perform clustering')
parser_clust.add_argument('input', type=str, help='input matrix in .npy format')
parser_clust.add_argument('output', type=str, help='output file')
parser_clust.add_argument('--coverage', help='coverage file')
parser_clust.set_defaults(func=sc.clustering)

# create the parser for the "embed" command
parser_embed = subparsers.add_parser('embed', help='Get UMAP embedding')
parser_embed.add_argument('input', type=str, help='input matrix in .npy format')
parser_embed.add_argument('output', type=str, help='output file')
parser_embed.set_defaults(func=sc.getEmbedding)

def main():
    args = parser.parse_args()
    args.func(args)