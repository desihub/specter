"""
2D Spectroperfectionism extractions in progress...
"""

import numpy as N
from scipy.sparse import spdiags, issparse

def sym_sqrt(a):
    """
    NAME: sym_sqrt

    PURPOSE: take 'square root' of a symmetric matrix via diagonalization

    USAGE: s = sym_sqrt(a)

    ARGUMENT: a: real symmetric square 2D ndarray

    RETURNS: s such that a = numpy.dot(s, s)

    WRITTEN: Adam S. Bolton, U. of Utah, 2009
    """
    
    w, v = N.linalg.eigh(a)
    w[w<0]=0 # Is this necessary to enforce eigenvalues positive definite???
        
    # dm = n.diagflat(n.sqrt(w))
    # result = n.dot(v, n.dot(dm, n.transpose(v)))

    #- A bit faster with sparse matrix for multiplication:
    nw = len(w)
    dm = spdiags(N.sqrt(w), 0, nw, nw)
    result = v.dot( dm.dot(v.T) )
    
    return result

def resolution_from_icov(icov):
    """
    Function to generate the 'resolution matrix' in the simplest
    (no unrelated crosstalk) Bolton & Schlegel 2010 sense.
    Works on dense matrices.  May not be suited for production-scale
    determination in a spectro extraction pipeline.

    Input argument is inverse covariance matrix array.
    If input is not 2D and symmetric, results will be unpredictable.

    WRITTEN: Adam S. Bolton, U. of Utah, 2009
    """
    if issparse(icov):
        icov = icov.toarray()
        
    sqrt_icov = sym_sqrt(icov)
    norm_vector = N.sum(sqrt_icov, axis=1)
    r_mat = N.outer(norm_vector**(-1), N.ones(norm_vector.size)) * sqrt_icov
    return r_mat
