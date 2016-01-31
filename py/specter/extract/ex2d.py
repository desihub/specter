"""
2D Spectroperfectionism extractions
"""

import sys
import numpy as np
import scipy.sparse
import scipy.linalg
from scipy.sparse import spdiags, issparse
from scipy.sparse.linalg import spsolve


def ex2d(image, ivar, psf, specrange, wavelengths, xyrange=None,
         full_output=False, regularize=0.0, ndecorr=False):
    """
    2D PSF extraction of flux from image given pixel inverse variance.
    
    Inputs:
        image : 2D array of pixels
        ivar  : 2D array of inverse variance for the image
        psf   : PSF object
        specrange : (start, stop) to extract
            (python style indexing; stop not included)
        wavelengths : 1D array of wavelengths to extract
        
    Optional Inputs:
        xyrange = (xmin, xmax, ymin, ymax): treat image as a subimage
            cutout of this region from the full image
        full_output : if True, return a dictionary of outputs including
            intermediate outputs such as the projection matrix.
        ndecorr : if True, decorrelate the noise between fibers, at the
            cost of residual signal correlations between fibers.
        
    Returns (flux, ivar, R):
        flux[nspec, nwave] = extracted resolution convolved flux
        ivar[nspec, nwave] = inverse variance of flux
        R : 2D resolution matrix to convert
    """

    #- Range of image to consider
    waverange = (wavelengths[0], wavelengths[-1])
    
    if xyrange is None:
        xmin, xmax, ymin, ymax = xyrange = psf.xyrange(specrange, waverange)
        image = image[ymin:ymax, xmin:xmax]
        ivar = ivar[ymin:ymax, xmin:xmax]
    else:
        xmin, xmax, ymin, ymax = xyrange

    nx, ny = xmax-xmin, ymax-ymin
    npix = nx*ny
    
    nspec = specrange[1] - specrange[0]
    nwave = len(wavelengths)
    
    #- Solve AT W pix = (AT W A) flux
    
    #- Projection matrix and inverse covariance
    A = psf.projection_matrix(specrange, wavelengths, xyrange)

    #- Pixel weights matrix
    w = ivar.ravel()
    W = spdiags(ivar.ravel(), 0, npix, npix)

    #-----
    #- Extend A with an optional regularization term to limit ringing.
    #- If any flux bins don't contribute to these pixels,
    #- also use this term to constrain those flux bins to 0.
    
    #- Original: exclude flux bins with 0 pixels contributing
    # ibad = (A.sum(axis=0).A == 0)[0]
    
    #- Identify fluxes with very low weights of pixels contributing            
    fluxweight = W.dot(A).sum(axis=0).A[0]
    minweight = 0.01*np.max(fluxweight)
    ibad = fluxweight < minweight
    
    #- Original version; doesn't work on older versions of scipy
    # I = regularize*scipy.sparse.identity(nspec*nwave)
    # I.data[0,ibad] = minweight - fluxweight[ibad]
    
    #- Add regularization of low weight fluxes
    Idiag = regularize*np.ones(nspec*nwave)
    Idiag[ibad] = minweight - fluxweight[ibad]
    I = scipy.sparse.identity(nspec*nwave)
    I.setdiag(Idiag)

    #- Only need to extend A if regularization is non-zero
    if np.any(I.diagonal()):
        pix = np.concatenate( (image.ravel(), np.zeros(nspec*nwave)) )
        Ax = scipy.sparse.vstack( (A, I) )
        wx = np.concatenate( (w, np.ones(nspec*nwave)) )
    else:
        pix = image.ravel()
        Ax = A
        wx = w

    #- Inverse covariance
    Wx = spdiags(wx, 0, len(wx), len(wx))
    iCov = Ax.T.dot(Wx.dot(Ax))

    #- Solve (image = A flux) weighted by Wx:
    #-     A^T W image = (A^T W A) flux = iCov flux    
    y = Ax.T.dot(Wx.dot(pix))
    
    xflux = spsolve(iCov, y).reshape((nspec, nwave))

    #- Solve for Resolution matrix
    try:
        if ndecorr:
            R, fluxivar = resolution_from_icov(iCov)
        else:
            R, fluxivar = resolution_from_icov(iCov, decorr=[nwave for x in range(nspec)])
    except np.linalg.linalg.LinAlgError, err:
        outfile = 'LinAlgError_{}-{}_{}-{}.fits'.format(specrange[0], specrange[1], waverange[0], waverange[1])
        print "ERROR: Linear Algebra didn't converge"
        print "Dumping {} for debugging".format(outfile)
        from astropy.io import fits
        fits.writeto(outfile, image, clobber=True)
        fits.append(outfile, ivar, name='IVAR')
        fits.append(outfile, A.data, name='ADATA') 
        fits.append(outfile, A.indices, name='AINDICES')
        fits.append(outfile, A.indptr, name='AINDPTR')
        fits.append(outfile, iCov.toarray(), name='ICOV')
        raise err
        
    #- Convolve with Resolution matrix to decorrelate errors
    fluxivar = fluxivar.reshape((nspec, nwave))
    rflux = R.dot(xflux.ravel()).reshape(xflux.shape)

    if full_output:
        results = dict(flux=rflux, ivar=fluxivar, R=R, xflux=xflux, A=A)
        results['iCov'] = iCov
        return results
    else:
        return rflux, fluxivar, R


def eigen_compose(w, v, invert=False, sqr=False):
    """
    Create a matrix from its eigenvectors and eigenvalues.

    Given the eigendecomposition of a matrix, recompose this
    into a real symmetric matrix.  Optionally take the square
    root of the eigenvalues and / or invert the eigenvalues.
    The eigenvalues are regularized such that the condition 
    number remains within machine precision for 64bit floating 
    point values.

    Args:
        w (array): 1D array of eigenvalues
        v (array): 2D array of eigenvectors.
        invert (bool): Should the eigenvalues be inverted? (False)
        sqr (bool): Should the square root eigenvalues be used? (False)

    Returns:
        A 2D numpy array which is the recomposed matrix.
    """
    dim = w.shape[0]

    # Threshold is 10 times the machine precision (~1e-15)
    threshold = 10.0 * sys.float_info.epsilon

    maxval = np.max(w)
    wscaled = np.zeros_like(w)

    if invert:
        # Normally, one should avoid explicit loops in python.
        # in this case however, we need to conditionally invert
        # the eigenvalues only if they are above the threshold.
        # Otherwise we might divide by zero.  Since the number
        # of eigenvalues is never too large, this should be fine.
        # If it does impact performance, we can improve this in
        # the future.  NOTE: simple timing with an average over
        # 10 loops shows that all four permutations of invert and
        # sqr options take about the same time- so this is not
        # an issue.
        if sqr:
            minval = np.sqrt(maxval) * threshold
            replace = 1.0 / minval
            tempsqr = np.sqrt(w)
            for i in range(dim):
                if tempsqr[i] > minval:
                    wscaled[i] = 1.0 / tempsqr[i]
                else:
                    wscaled[i] = replace
        else:
            minval = maxval * threshold
            replace = 1.0 / minval
            for i in range(dim):
                if w[i] > minval:
                    wscaled[i] = 1.0 / w[i]
                else:
                    wscaled[i] = replace
    else:
        if sqr:
            minval = np.sqrt(maxval) * threshold
            replace = minval
            wscaled[:] = np.where((w > minval), np.sqrt(w), replace*np.ones_like(w))
        else:
            minval = maxval * threshold
            replace = minval
            wscaled[:] = np.where((w > minval), w, replace*np.ones_like(w))

    # multiply to get result
    wdiag = spdiags(wscaled, 0, dim, dim)
    return v.dot( wdiag.dot(v.T) )


def resolution_from_icov(icov, decorr=None):
    """
    Function to generate the 'resolution matrix' in the simplest
    (no unrelated crosstalk) Bolton & Schlegel 2010 sense.
    Works on dense matrices.  May not be suited for production-scale
    determination in a spectro extraction pipeline.

    Args:
        icov (array): real, symmetric, 2D array containing inverse
                      covariance.
        decorr (list): produce a resolution matrix which decorrelates
                      signal between fibers, at the cost of correlated
                      noise between fibers (default).  This list should
                      contain the number of elements in each spectrum,
                      which is used to define the size of the blocks.

    Returns (R, ivar):
        R : resolution matrix
        ivar : R C R.T  -- decorrelated resolution convolved inverse variance
    """
    #- force symmetry since due to rounding it might not be exactly symmetric
    icov = 0.5*(icov + icov.T)
    
    if issparse(icov):
        icov = icov.toarray()

    w, v = scipy.linalg.eigh(icov)

    sqrt_icov = np.zeros_like(icov)

    if decorr is not None:
        if np.sum(decorr) != icov.shape[0]:
            raise RuntimeError("The list of spectral block sizes must sum to the matrix size")
        inverse = eigen_compose(w, v, invert=True)
        # take each spectrum block and process
        offset = 0
        for b in decorr:
            bw, bv = scipy.linalg.eigh(inverse[offset:offset+b,offset:offset+b])
            sqrt_icov[offset:offset+b,offset:offset+b] = eigen_compose(bw, bv, invert=True, sqr=True)
            offset += b
    else:
        sqrt_icov = eigen_compose(w, v, sqr=True)

    norm_vector = np.sum(sqrt_icov, axis=1)
    R = np.outer(norm_vector**(-1), np.ones(norm_vector.size)) * sqrt_icov
    ivar = norm_vector**2  #- Bolton & Schlegel 2010 Eqn 13
    return R, ivar
