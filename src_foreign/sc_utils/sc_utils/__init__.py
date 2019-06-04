import scipy as sp
import numpy as np
import gzip
import math

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

'''
def ldaTransform(args):
    data = InputData(args.input)
    n_dim = 30

    lda = LdaMulticore(data, num_topics=n_dim)
'''


def lsiTransform(args):
    from gensim.models.lsimodel import stochastic_svd
    from gensim.matutils import corpus2dense, corpus2csc
    from gensim.models import LsiModel

    data = InputData(args.input)
    n_dim = 30

    model = LsiModel(data, num_topics=n_dim, chunksize=2000)
    data_transformed = corpus2dense(model[data], n_dim).T
    np.savetxt(args.output, data_transformed, delimiter='\t')

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
    import umap

    mat = np.loadtxt(args.input)
    embedding = umap.UMAP(random_state=42, n_components=3).fit_transform(mat)
    np.savetxt(args.output, embedding, delimiter='\t')

# regress out a variable
def regressOut(X, y):
    from sklearn.linear_model import Lasso
    reg = Lasso(alpha=0.1).fit(X, y)
    return y - reg.predict(X)

def clustering(args):
    from sklearn.neighbors import kneighbors_graph
    import igraph as ig
    import leidenalg as la

    def readCoverage(fl):
        with open(fl, 'r') as f:
            return [[math.log(int(line.strip()))] for line in f]

    data_transformed = np.loadtxt(args.input)
    print(data_transformed.shape)
    if (args.coverage):
        print("Performing regression")
        cov = readCoverage(args.coverage)
        def normalize(y):
            return [r[0] for r in regressOut(cov, [[x] for x in y])]
        data_transformed = np.apply_along_axis(normalize, 0, data_transformed)
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