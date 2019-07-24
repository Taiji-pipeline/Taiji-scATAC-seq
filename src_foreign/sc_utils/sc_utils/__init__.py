import scipy as sp
import numpy as np
import math
from gensim.matutils import corpus2dense, corpus2csc

from .LSI import lsiTransform
from .DiffusionMap import diffusionMap

'''
    def transform(self, u):
        for doc in self:
            vec = corpus2csc([doc], num_terms=self.num_terms)
            yield vec.T * u
'''

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
    if (args.dim):
        data_transformed = data_transformed[..., :args.dim]
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
    adj = kneighbors_graph(data_transformed, args.k, mode='distance')

    print("Start clustering")
    vcount = max(adj.shape)
    sources, targets = adj.nonzero()
    edgelist = list(zip(sources.tolist(), targets.tolist()))
    gr = ig.Graph(vcount, edgelist)

    if(args.res):
        partition = la.find_partition(gr, la.CPMVertexPartition,
            n_iterations=10, seed=12343, resolution_parameter = args.res)
    else:
        partition = la.find_partition(gr, la.ModularityVertexPartition,
            n_iterations=10, seed=12343)

    print("Clusters: ")
    print(len(partition)) 
    with open(args.output, 'w') as f:
        f.write('\n'.join([','.join(map(str, c)) for c in partition]))

    if(args.embed):
        print("Create Embedding:")
        getEmbedding(data_transformed, args.embed, method=args.embed_method)
