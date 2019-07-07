import scipy as sp
import numpy as np
import gzip
import math
from gensim.matutils import corpus2dense, corpus2csc

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

def diffusionMap(args):
    data = InputData(args.input)
    n_dim = 15
    t = 1
    print("Read Data")
    indptr = [0]
    indices = []
    mat = []
    for row in iter(data):
        for (i,x) in row:
            indices.append(i)
            mat.append(1)
        indptr.append(len(indices))
    mat = sp.sparse.csr_matrix((mat, indices, indptr), dtype=int)
    (n,_) = mat.get_shape()

    print("Compute similarity matrix")
    jm = mat.dot(mat.T).todense()
    s = mat.sum(axis=1).dot(np.ones((1,n)))
    jm = jm / (s + s.T - jm)
    np.fill_diagonal(jm, 0)
    print("Normalization")
    s = jm.sum(axis=1).dot(np.ones((1,n)))
    jm = jm / s

    print("Reduction")
    (evals, evecs) = sp.sparse.linalg.eigs(jm, k=n_dim+1, which='LR')
    ix = evals.argsort()[::-1][1:]
    evals = np.real(evals[ix])
    evecs = np.real(evecs[:, ix])
    dmap = np.matmul(evecs, np.diag(evals**t))
    np.savetxt(args.output, dmap, delimiter='\t')


def ldaTransform(args):
    from gensim.models.ldamodel import LdaModel
    import logging
    logging.basicConfig(filename='gensim.log',
        format="%(asctime)s:%(levelname)s:%(message)s",
        level=logging.INFO)

    data = InputData(args.input)
    n_dim = 15

    model = LdaModel(data, num_topics=n_dim, chunksize=10000, random_state=2347, passes=40, iterations=5000, eval_every=10)
    data_transformed = corpus2dense(model[data], n_dim).T
    np.savetxt(args.output, data_transformed, delimiter='\t')

def lsiTransform(args):
    from gensim.models.lsimodel import stochastic_svd
    from gensim.models import LsiModel

    data = InputData(args.input)
    n_dim = 15

    model = LsiModel(data, num_topics=n_dim, onepass=False, power_iters=2, chunksize=10000)
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

def reduceDimension(args):
    if(args.method == "svd"):
        lsiTransform(args)
    elif(args.method == "dm"):
        diffusionMap(args)
    else:
        ldaTransform(args)

def getEmbedding(mat, output, method="umap"):
    print(method)
    if(method == "umap"):
        import umap
        e1 = umap.UMAP(random_state=42,
            n_components=2, min_dist=0).fit_transform(mat)
        e2 = umap.UMAP(random_state=42,
            n_components=3, min_dist=0).fit_transform(mat)
        embedding = np.concatenate((e1, e2), axis=1)
    elif(method == "none"):
        e1 = mat[...,:2]
        e2 = mat[...,:3]
        embedding = np.concatenate((e1, e2), axis=1)
    else:
        from MulticoreTSNE import MulticoreTSNE as TSNE
        e1 = TSNE(n_jobs=4, n_components=2, perplexity=30).fit_transform(mat)
        e2 = TSNE(n_jobs=4, n_components=3, perplexity=30).fit_transform(mat)
        embedding = np.concatenate((e1, e2), axis=1)
    np.savetxt(output, embedding, delimiter='\t')

# regress out a variable
def regressOut(X, y):
    from sklearn.linear_model import Lasso
    reg = Lasso(alpha=0.1).fit(X, y)
    return np.subtract(y.flatten(), reg.predict(X))

def clustering(args):
    from sklearn.neighbors import kneighbors_graph
    import igraph as ig
    import leidenalg as la

    def readCoverage(fl):
        with open(fl, 'r') as f:
            return np.array([[math.log(int(line.strip()))] for line in f])

    def scaling(xs):
        s = 0
        for x in xs:
            s = s + x * x
        s = math.sqrt(s)
        return np.array([x / s for x in xs])

    data_transformed = np.loadtxt(args.input)
    print(data_transformed.shape)

    if (args.coverage):
        print("Performing regression")
        cov = readCoverage(args.coverage)
        def normalize(y):
            return regressOut(cov, np.array([[x] for x in y]))
        data_transformed = np.apply_along_axis(normalize, 0, data_transformed)

    if (args.discard):
        data_transformed = data_transformed[..., 1:]

    if (args.scale):
        data_transformed = np.apply_along_axis(scaling, 1, data_transformed)

    print(data_transformed.shape)

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

    if(args.embed):
        print("Create Embedding:")
        getEmbedding(data_transformed, args.embed, method=args.embed_method)
