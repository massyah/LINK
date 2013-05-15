
import logging
import math
import os
import itertools

import numpy
from numpy import dot as numdot
from scipy import dot as scidot
import scipy.sparse
import scipy.linalg
from scipy.linalg.lapack import get_lapack_funcs, find_best_lapack_type
from gensim import *
from gensim.matutils import cossim

# get_ipython().magic("cd ../material")
# get_ipython().magic("run -i svd_factorization.py")
# load_corpus("33902")
# get_ipython().magic("cd ../SACE\ analysis/")


# scipy is not a stable package yet, locations change, so try to work
# around differences (currently only concerns location of 'triu' in scipy 0.7 vs. 0.8)
try:
    from scipy.linalg.basic import triu
except ImportError:
    from scipy.linalg.special_matrices import triu

blas = lambda name, ndarray: scipy.linalg.get_blas_funcs((name,), (ndarray,))[0]

blas_nrm2 = blas('nrm2', numpy.array([], dtype = float))
blas_scal = blas('scal', numpy.array([], dtype = float))





def cus_unitVec(vec):
    """
    Scale a vector to unit length. The only exception is the zero vector, which
    is returned back unchanged.
    
    If the input is sparse (list of 2-tuples), output will also be sparse. Otherwise,
    output will be a numpy array.
    """
    if scipy.sparse.issparse(vec): # convert scipy.sparse to standard numpy array
        vec = vec.toarray().flatten()
    
    try:
        first = iter(vec).next() # is there at least one element?
    except:
        return vec
    
    if isinstance(first, tuple): # sparse format?
        vecLen = 1.0 * math.sqrt(sum(val * val for _, val in vec))
        assert vecLen > 0.0, "sparse documents must not contain any explicit zero entries"
        if vecLen != 1.0:
            return [(termId, val / vecLen) for termId, val in vec]
        else:
            return list(vec)
    else: # dense format
        vec = numpy.asarray(vec, dtype=float)
        veclen = blas_nrm2(vec)
        if veclen > 0.0:
            return blas_scal(1.0 / veclen, vec)
        else:
            return vec


def cus_cossim(vec1, vec2):
    # vec1, vec2 = dict(vec1), dict(vec2)
    if not vec1 or not vec2:
        return 0.0
    # vec1Len = 1.0 * math.sqrt(sum(val * val for val in vec1.itervalues()))
    # vec2Len = 1.0 * math.sqrt(sum(val * val for val in vec2.itervalues()))
    # assert vec1Len > 0.0 and vec2Len > 0.0, "sparse documents must not contain any explicit zero entries"
    if len(vec2) < len(vec1):
        vec1, vec2 = vec2, vec1 # swap references so that we iterate over the shorter vector
    result = sum(value * vec2.get(index, 0.0) for index, value in vec1.iteritems())
    # result /= vec1Len * vec2Len # rescale by vector lengths
    return result


NTERMS=7967821



vec1=[(1089, 1), (2002, 1), (2145, 26), (2166, 5), (2341, 1), (3347, 1), (3384, 1), (4075, 1), (5569, 2), (6051, 1), (6081, 1), (6713, 1), (7878, 3), (7958, 1), (8424, 1), (10805, 1), (11174, 1), (11242, 4), (11819, 2), (12266, 1), (14266, 2), (14581, 1), (15155, 2), (16371, 1), (17365, 1), (17419, 1), (17481, 1), (17737, 1), (19065, 5), (19447, 1), (19861, 1), (20113, 1), (20114, 1), (20175, 1), (20356, 1), (20721, 1), (20726, 1), (20828, 1), (20924, 1), (21428, 1), (21625, 2), (22458, 1), (22942, 1), (23006, 1), (23633, 1), (24041, 1), (24263, 1), (26368, 1), (26424, 1), (26777, 1), (27408, 1), (27481, 2), (27890, 2), (27892, 1), (28230, 3), (28621, 1), (32071, 1), (32469, 1), (33104, 1), (33484, 1), (33504, 2), (34453, 2), (34804, 1), (36394, 1), (36975, 1), (37059, 2), (37345, 1), (38890, 1), (39081, 4), (40066, 1), (41857, 1), (41982, 1), (43110, 1), (43450, 1), (43692, 1), (43923, 1), (44001, 2), (44154, 1), (44213, 1), (44883, 1), (45627, 1), (45738, 2), (46736, 2), (47340, 2), (47901, 1), (48515, 1), (49210, 1), (50558, 1), (50560, 1), (50697, 1), (51477, 1), (52271, 1), (52845, 1), (53677, 2), (56245, 1), (56556, 1), (56808, 1), (57443, 4), (57539, 1), (57556, 2), (58288, 2), (58783, 4), (60439, 1), (61025, 1), (62243, 2), (63997, 1), (64407, 4), (66864, 1), (67249, 2), (68036, 1), (68049, 1), (70529, 2), (72215, 1), (72915, 1), (73804, 1)]
vec2=[(586, 2), (1376, 2), (1598, 1), (2166, 2), (3285, 1), (5096, 1), (6081, 1), (9971, 2), (11242, 2), (14266, 1), (14581, 1), (15783, 1), (17060, 1), (17491, 1), (18381, 2), (19244, 4), (20356, 1), (20721, 2), (20726, 2), (21625, 1), (21854, 1), (24287, 2), (26424, 2), (27889, 1), (28230, 1), (30371, 1), (31078, 1), (31985, 2), (32071, 1), (33416, 1), (34281, 1), (34804, 1), (35634, 1), (36662, 2), (36968, 1), (39226, 4), (40214, 1), (41696, 2), (41706, 1), (42065, 6), (42358, 1), (43110, 1), (43923, 1), (44001, 1), (46538, 1), (48473, 1), (48515, 1), (49067, 1), (50697, 1), (51280, 1), (51911, 2), (53677, 2), (56245, 1), (56323, 1), (56808, 1), (57443, 1), (57579, 1), (58288, 1), (58783, 1), (59727, 1), (61037, 14), (63454, 6), (63710, 2), (64407, 5), (64820, 1), (67962, 1), (68036, 1), (68413, 1), (68454, 2), (69330, 2), (70829, 1), (72097, 1), (72501, 1), (72789, 2)]

vec1,vec2=tfidf[vec1],tfidf[vec2]

vec1lsi=sparse2full(cus_unitVec(lsi[tfidf[vec1]]),500)
vec2lsi=sparse2full(cus_unitVec(lsi[tfidf[vec2]]),500)


## We have the equivalence between
# cus_unitVec(lsi[tfidf[corpus[0]]])
# #and
# cus_unitvec(index.corpus[0])


print cossim(vec1,vec2)
#Possible opt to cus_cossim: 
## Use normalized vec as input
vec1,vec2=cus_unitVec(vec1),cus_unitVec(vec2)
## Use dict as input
vec1,vec2=dict(vec1),dict(vec2)
## 
assert(abs(cus_cossim(vec1,vec2)-0.0112255821) < 0.001)
# Build the scipy sparse rep

vec1sp=scipy.sparse.dok_matrix((1,NTERMS),dtype=float)
for k,v in vec1.items():
	vec1sp[0,k]=v

vec2sp=scipy.sparse.dok_matrix((1,NTERMS),dtype=float)
for k,v in vec2.items():
	vec2sp[0,k]=v

vec1dense=vec1sp.todense()
vec2dense=vec2sp.todense()

get_ipython().magic("timeit numdot(vec1dense,vec1dense.T)[0,0]")
get_ipython().magic("timeit numdot(vec1lsi,vec2lsi.T)")
get_ipython().magic("timeit cus_cossim(vec1,vec1)")
get_ipython().magic("timeit cossim(vec1,vec1)")
# get_ipython().magic("timeit numdot(vec1sp,vec1sp.T)[0,0]")
