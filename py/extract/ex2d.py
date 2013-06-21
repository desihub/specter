"""
2D Spectroperfectionism extractions in progress...
"""

import numpy as N
from scipy.sparse import spdiags, issparse
from scipy.sparse.linalg import spsolve

def ex2d(image, ivar, psf, specrange, wavelengths, full_output=False):
    """
    2D PSF extraction of flux from image given pixel inverse variance.
    
    Inputs:
        image : 2D array of pixels
        ivar  : 2D array of inverse variance for the image
        psf   : PSF object
        specrange : (specmin, specmax) inclusive to extract
        wavelengths : 1D array of wavelengths to extract
        
    Returns (flux, ivar, R):
        flux[nspec, nwave] = extracted resolution convolved flux
        ivar[nspec, nwave] = inverse variance of flux
        R : 2D resolution matrix to convert
        
    if full_output is True, return dictionary of the above plus additional
    intermediate outputs such as the projection matrix.
        
    To do: Add ability to input subimage instead of entire image
    """
    
    #- Range of image to consider
    waverange = (wavelengths[0], wavelengths[-1])
    xmin, xmax, ymin, ymax = xyrange = psf.xyrange(specrange, waverange)
    nx, ny = xmax-xmin, ymax-ymin
    npix = nx*ny
    image = image[ymin:ymax, xmin:xmax]
    ivar = ivar[ymin:ymax, xmin:xmax]
    
    nspec = specrange[1] - specrange[0]
    nflux = len(wavelengths)
    
    #- Pixel weights matrix
    W = spdiags(ivar.ravel(), 0, npix, npix)
    
    #- Projection matrix and inverse covariance
    A = psf.projection_matrix(specrange, wavelengths, xyrange)
    iCov = A.T.dot(W.dot(A))
    
    #- Solve image = A flux weighted by W:
    #-     A^T W image = (A^T W A) flux = iCov flux    
    y = A.T.dot(W.dot(image.ravel()))
    
    xflux = spsolve(iCov, y).reshape((nspec, nflux))

    #- Convolve with Resolution matrix to decorrelate errors
    R = resolution_from_icov(iCov)
    rflux = R.dot(xflux.ravel()).reshape(xflux.shape)

    #- Note: this could be faster by using calculations already done in
    #- resolution_from_iCov()
    Cov = N.linalg.inv(iCov.toarray())
    rCov = R.T.dot(Cov.dot(R))
    ivar = 1.0/rCov.diagonal().reshape(rflux.shape)

    if full_output:
        results = dict(flux=rflux, ivar=ivar, R=R, xflux=xflux, projmat=A)
        results['Cov'] = Cov
        results['iCov'] = iCov
        results['rCov'] = rCov
        return results
    else:
        return rflux, ivar, R
    

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
