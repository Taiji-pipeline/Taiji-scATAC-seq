import scipy as sp
import numpy as np
from gensim.models.lsimodel import stochastic_svd
from gensim.models import LsiModel
from gensim.matutils import corpus2dense, corpus2csc

from .Utils import InputData

def lsiTransform(args):

    data = InputData(args.input)
    n_dim = 15

    #model = LsiModel(data, num_topics=n_dim, onepass=False, power_iters=2, chunksize=10000)
    model = LsiModel(data, num_topics=n_dim, chunksize=10000)
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

