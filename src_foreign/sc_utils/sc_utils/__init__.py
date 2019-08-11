import scipy as sp
import numpy as np
import math
from gensim.matutils import corpus2dense, corpus2csc
from sklearn.neighbors import kneighbors_graph
import igraph as ig

from .DiffusionMap import diffusionMap

def reduceDimension(args):
    if(args.method == "svd"):
        lsiTransform(args)
    else:
        diffusionMap(args)

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

def clustering(args):
    import leidenalg as la

    data = readCoordinates(args.input, n_dim=args.dim,
        discard=args.discard, scale=args.scale)

    print("Start KNN")
    gr = mkKNNGraph([args.input], k=args.k)

    print("Start clustering")
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
        getEmbedding(data, args.embed, method=args.embed_method)

def readCoordinates(fl, n_dim=None, discard=None, scale=None):
    def scaling(xs):
        s = 0
        for x in xs:
            s = s + x * x
        s = math.sqrt(s)
        return np.array([x / s for x in xs])
    data = np.loadtxt(fl)
    if (n_dim):
        data = data[..., :n_dim]
    if (discard):
        data = data[..., 1:]
    if (scale):
        data = np.apply_along_axis(scaling, 1, data)
    return data

def mkKNNGraph(fls, k=25):
    adj = None
    for fl in fls:
        mat = readCoordinates(fl)
        if (adj == None):
            adj = kneighbors_graph(mat, k, mode='distance')
        else:
            adj += kneighbors_graph(mat, k, mode='distance')
    vcount = max(adj.shape)
    sources, targets = adj.nonzero()
    edgelist = list(zip(sources.tolist(), targets.tolist()))
    gr = ig.Graph(vcount, edgelist)
    return gr