import argparse
import scipy as sp
import numpy as np
import gzip
import umap
from gensim.models.lsimodel import stochastic_svd
from gensim.matutils import corpus2dense, corpus2csc
from gensim.models import LsiModel
from sklearn.neighbors import kneighbors_graph
import igraph as ig
import leidenalg as la

from sklearn.utils.extmath import randomized_svd

class InputData:
    def __init__(self, filename):
        self.filename = filename
        with gzip.open(self.filename, mode='rt') as f:
            header = f.readline().strip()
            (m, n) = header.split(": ")[1].split(" x ")
            self.num_doc = int(m)
            self.num_terms = int(n)

    def __iter__(self):
        def convert(x):
            return (int(x[0]), float(x[1]))
        with gzip.open(self.filename, mode='rt') as f:
            f.readline()
            for line in f:
                yield [convert(item.split(",")) for item in line.strip().split("\t")[1:]]

'''
    def transform(self, u):
        for doc in self:
            vec = corpus2csc([doc], num_terms=self.num_terms)
            yield vec.T * u
'''

def lsiTransform(args):
    data = InputData(args.input)
    n_dim = 30

    model = LsiModel(data, num_topics=n_dim, chunksize=2000)
    data_transformed = corpus2dense(model[data], n_dim).T
    np.save(args.output, data_transformed)

    ''' SVD
    (u, s) = stochastic_svd(data, n_dim, data.num_terms, power_iters=2, chunksize=2000)
    data_transformed = sp.concatenate(list(data.transform(u)))
    '''

    ''' scipy SVD
    data2 = corpus2dense(data, data.num_terms)
    (u2, s2, v2) = randomized_svd(data2, 
                                n_components=5,
                                n_iter=5,
                                random_state=None)
    print((sp.sparse.diags(s2) * v2).T)
    '''


def getEmbedding(args):
    mat = np.load(args.input)
    embedding = umap.UMAP(random_state=42).fit_transform(mat)
    np.savetxt(args.output, embedding, delimiter='\t')

def clustering(args):
    data_transformed = np.load(args.input)
    print("Start KNN")
    adj = kneighbors_graph(data_transformed, 20, mode='distance')

    print("Start clustering")
    vcount = max(adj.shape)
    sources, targets = adj.nonzero()
    edgelist = list(zip(sources.tolist(), targets.tolist()))
    gr = ig.Graph(vcount, edgelist)

    partition = la.find_partition(gr, la.ModularityVertexPartition, n_iterations=-1, seed=12343)
    print("Clusters: ")
    print(len(partition)) 
    with open(args.output, 'w') as f:
        f.write('\n'.join([','.join(map(str, c)) for c in partition]))

################################################################################
## ARGUMENT PARSER
################################################################################

parser = argparse.ArgumentParser(description="Fast online SVD")
subparsers = parser.add_subparsers(title="sub commands")

# create the parser for the "run" command
parser_run = subparsers.add_parser('run', help='run SVD')
parser_run.add_argument('input', type=str, help='gzipped input file')
parser_run.add_argument('output', type=str, help='output matrix in .npy format')
parser_run.set_defaults(func=lsiTransform)

# create the parser for the "clust" command
parser_clust = subparsers.add_parser('clust', help='perform clustering')
parser_clust.add_argument('input', type=str, help='input matrix in .npy format')
parser_clust.add_argument('output', type=str, help='output file')
parser_clust.set_defaults(func=clustering)

# create the parser for the "embed" command
parser_embed = subparsers.add_parser('embed', help='Get UMAP embedding')
parser_embed.add_argument('input', type=str, help='input matrix in .npy format')
parser_embed.add_argument('output', type=str, help='output file')
parser_embed.set_defaults(func=getEmbedding)

args = parser.parse_args()
args.func(args)
