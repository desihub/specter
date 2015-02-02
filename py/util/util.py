"""
Utility functions and classes for specter

Stephen Bailey
Fall 2012
"""

import math
import numpy as N
import scipy.signal
from scipy.special import legendre
from scipy.sparse import spdiags
from scipy.signal import convolve, convolve2d
from specter.util import pixspline
from time import time

from specter.extract.ex2d import resolution_from_icov

_t0 = 0.0
def _timeit():
    global _t0
    tx = time()
    dt = tx - _t0
    _t0 = tx
    return dt
    
#- 2D Linear interpolator
class LinearInterp2D(object):
    """
    Linear interpolation on a 2D grid.  Allows values to be interpolated
    to be multi-dimensional.
    """
    def __init__(self, x, y, data):
        """
        x : array of x coordinates
        y : array of y coordinates
        data[ix, iy, ...] : 3 or more dimensional array of data to interpolate
            first two coordinates are x and y
        """
        self.x = N.array(x)
        self.y = N.array(y)
        self.data = N.array(data)

    def __call__(self, x, y):
        """
        Evaluate data at (x,y)
        """
        #- TODO: compare speed to solution at
        #- http://stackoverflow.com/questions/12729228/simple-efficient-bilinear-interpolation-of-images-in-numpy-and-python
        
        #- Find where we are in grid
        #- clip to 1 because we will use i and i-1
        #- clip to len(x)-1 to allow extrapolation beyond grid boundary
        ix = N.searchsorted(self.x, x).clip(1, len(self.x)-1)
        iy = N.searchsorted(self.y, y).clip(1, len(self.y)-1)
        
        #- Interpolation distances from points
        dx = (x - self.x[ix-1]) / (self.x[ix] - self.x[ix-1])
        dy = (y - self.y[iy-1]) / (self.y[iy] - self.y[iy-1])

        #- Interpolate, allowing x and/or y to be multi-dimensional
        #- NOTE: these are the slow steps, about equal time each
        
        #- Original code with what appears to be vestigial transposes
        # data1 = (self.data[ix-1,iy-1].T*(1-dx) + self.data[ix,iy-1].T*dx).T
        # data2 = (self.data[ix-1,iy].T*(1-dx) + self.data[ix,iy].T*dx).T
        # dataxy = (data1.T*(1-dy) + data2.T*dy).T

        #- Updated without transposes
        data1 = (self.data[ix-1,iy-1]*(1-dx) + self.data[ix,iy-1]*dx)
        data2 = (self.data[ix-1,iy]*(1-dx) + self.data[ix,iy]*dx)
        dataxy = (data1*(1-dy) + data2*dy)

        return dataxy
        
def psfbias(p1, p2, wave, phot, ispec=0, readnoise=3.0):
    """
    Return bias from extracting with PSF p2 if the real PSF is p1
    
    Inputs:
        p1, p2 : PSF objects
        wave[] : wavelengths in Angstroms
        phot[] : spectrum in photons
        
    Optional Inputs:
        ispec : spectrum number
        readnoise : CCD read out noise (optional)
        
    Returns:
        bias array same length as wave
    """
    #- flux -> pixels projection matrices
    xyrange = p1.xyrange( (ispec,ispec+1), (wave[0], wave[-1]) )
    A = p1.projection_matrix((ispec,ispec+1), wave, xyrange)
    B = p2.projection_matrix((ispec,ispec+1), wave, xyrange)

    #- Pixel weights from photon shot noise and CCD read noise
    img = A.dot(phot)            #- True noiseless image
    imgvar = readnoise**2 + img  #- pixel variance
    npix = img.size
    W = spdiags(1.0/imgvar, 0, npix, npix)
    
    #- covariance matrix for each PSF
    iACov = A.T.dot(W.dot(A))
    iBCov = B.T.dot(W.dot(B))
    BCov = N.linalg.inv(iBCov.toarray())
    
    #- Resolution matricies
    RA, _ = resolution_from_icov(iACov)
    RB, _ = resolution_from_icov(iBCov)

    #- Bias
    bias = (RB.dot(BCov.dot(B.T.dot(W.dot(A)).toarray())) - RA).dot(phot) / RA.dot(phot)

    return bias

def psfabsbias(p1, p2, wave, phot, ispec=0, readnoise=3.0):
    """
    Return absolute bias from extracting with PSF p2 if the real PSF is p1.
    
    Inputs:
        p1, p2 : PSF objects
        wave[] : wavelengths in Angstroms
        phot[] : spectrum in photons
        
    Optional Inputs:
        ispec : spectrum number
        readnoise : CCD read out noise (optional)
        
    Returns bias, R
        bias array same length as wave
        R resolution matrix for PSF p1
        
        
    See psfbias() for relative bias
    """
    #- flux -> pixels projection matrices
    xyrange = p1.xyrange( (ispec,ispec+1), (wave[0], wave[-1]) )
    A = p1.projection_matrix((ispec,ispec+1), wave, xyrange)
    B = p2.projection_matrix((ispec,ispec+1), wave, xyrange)

    #- Pixel weights from photon shot noise and CCD read noise
    img = A.dot(phot)            #- True noiseless image
    imgvar = readnoise**2 + img  #- pixel variance
    npix = img.size
    W = spdiags(1.0/imgvar, 0, npix, npix)
    
    #- covariance matrix for each PSF
    iACov = A.T.dot(W.dot(A))
    iBCov = B.T.dot(W.dot(B))
    BCov = N.linalg.inv(iBCov.toarray())
    
    #- Resolution matricies
    RA, _ = resolution_from_icov(iACov)
    RB, _ = resolution_from_icov(iBCov)

    #- Bias
    bias = (RB.dot(BCov.dot(B.T.dot(W.dot(A)).toarray())) - RA).dot(phot)

    return bias, RA

def rebin_image(image, n):
    """
    rebin 2D array pix into bins of size n x n
    
    New binsize must be evenly divisible into original pix image
    """
    assert image.shape[0] % n == 0
    assert image.shape[1] % n == 0
    
    s = image.shape[0]//n, n, image.shape[1]//n, n
    return image.reshape(s).sum(-1).sum(1)

    
#- Utility functions for sinc shifting pixelated PSFs
def _sincfunc(x, dx, dampfac=3.25):
    """sinc helper function for sincshift()"""
    if dx != 0.0:
        xx = (x+dx)*N.pi  #- cache shifted array for 30% faster evals
        return N.exp( -(xx/(dampfac*N.pi))**2 ) * N.sin(xx) / xx
    else:
        xx = N.zeros(len(x))
        xx[len(x)/2] = 1.0
        return xx

#- Implementation note: the typical PSF image is 15x15.
#- fftconvolve is not faster than convolve for images this small
def sincshift(image, dx, dy, sincrad=10, dampfac=3.25):
    """
    Return image shifted by dx, dy using sinc interpolation.
    
    For speed, do each dimension independently which can introduce edge
    effects.  Also see sincshift2d().
    """
    s = N.arange(-sincrad, sincrad+1.0)
    imgshape = image.shape

    if abs(dx) > 1e-6:
        sincx = _sincfunc(s, -dx, dampfac=dampfac)
        image = convolve(image.ravel(), sincx, mode='same')
        image = image.reshape(imgshape)

    if abs(dy) > 1e-6:
        sincy = _sincfunc(s, -dy, dampfac=dampfac)
        image = convolve(image.T.ravel(), sincy, mode='same')
        image = image.reshape(imgshape[-1::-1]).T

    return image

def sincshift2d(image, dx, dy, sincrad=10, dampfac=3.25):
    """
    Return image shifted by dx, dy using full 2D sinc interpolation
    """
    s = N.arange(-sincrad, sincrad+1.0)
    sincx = _sincfunc(s, -dx, dampfac=dampfac)
    sincy = _sincfunc(s, -dy, dampfac=dampfac)
    kernel = N.outer(sincy, sincx)
    newimage = convolve2d(image, kernel, mode='same')
    return newimage

from scipy.special import erf
def gaussint(x, mean=0.0, sigma=1.0):
    """
    Return integral from -inf to x of normalized Gaussian with mean and sigma
    """
    z = (x - mean) / (math.sqrt(2) * sigma)
    return (erf(z) + 1.0) / 2.0
    
def gausspix(x, mean=0.0, sigma=1.0):
    """
    Return Gaussian(mean,sigma) integrated over unit width pixels centered at x[].
    """
    edges = N.concatenate((x-0.5, x[-1:]+0.5))
    integrals = gaussint(edges, mean=mean, sigma=sigma)
    return integrals[1:] - integrals[0:-1]
    
def weighted_solve(A, b, w):
    """
    Solve `A x = b` with weights `w` on `b`
    Returns x, inverseCovarance(x)
    """
    assert len(b) == len(w)
    
    n = len(b)
    W = spdiags(w, [0,], n, n)
    y = A.T.dot(W.dot(b))
    iCov = A.T.dot(W.dot(A))
    x = N.linalg.lstsq(iCov, y)[0]
    return x, iCov

def trapz(edges, xp, yp):
    """
    Perform trapezoidal integration between edges using sampled function
    yp vs. xp.  Returns array of length len(edges)-1.
        
    Input xp array must be sorted in ascending order.
    
    See also numpy.trapz, which integrates a single array
    """
    if N.any(N.diff(xp) < 0.0):
        raise ValueError("Input x must be sorted in increasing order")
    
    if len(xp) != len(yp):
        raise ValueError("xp and yp must have same length")

    yedge = N.interp(edges, xp, yp)
    result = N.zeros(len(edges)-1)
    iedge = N.searchsorted(xp, edges)
    for i in range(len(edges)-1):
        ilo, ihi = iedge[i], iedge[i+1]
        xx = N.concatenate( (edges[i:i+1], xp[ilo:ihi], edges[i+1:i+2]) )
        yy = N.concatenate( (yedge[i:i+1], yp[ilo:ihi], yedge[i+1:i+2]) )
        result[i] = N.trapz(yy, xx)
        
    return result
    
def model_function(x, y):
    """
    Find the function p(x) which integrates to y over each bin.
    
    Current implementation uses linear trapz interpolation between points.
    
    TODO: Refactor this name and interface.  Document underlying math.
    TODO: Replace with quadratic spline instead of linear.
    """
    nx = len(x)
    B = N.ones((3, nx)) * 0.75  #- B[1] = 0.75 for all entries
    B[0,1] = B[2, -2] = 0.25
    B[0,0] = B[-1,-1] = 0.0     #- doesn't really matter; unused
    dx1 = x[1:] - x[0:-1]
    dx2 = x[2:] - x[0:-2]
    # for i in range(1, nx-1):
    #     B[0,i+1] = 0.25*dx1[i]/dx2[i-1]
    #     B[2,i-1] = 0.25*dx1[i-1]/dx2[i-1]    
    B[0,2:nx] = 0.25 * dx1[1:nx-1] / dx2[0:nx-2]
    B[2,0:nx-2] = 0.25 * dx1[0:nx-2] / dx2[0:nx-2]
        
    p = scipy.linalg.solve_banded((1,1), B, y)
    return p

def get_bin_edges(bin_centers):
    """
    Given a sorted array of bin centers, return array of bin edges
    
    edge[0] = center[0]
    edge[-1] = center[-1]
    edge[i] = 0.5*(center[i-1] + center[i]) for 0 < i < len(center)-1
    """
    mid = 0.5*(bin_centers[0:-1] + bin_centers[1:])
    return N.concatenate( (bin_centers[0:1], mid, bin_centers[-1:]) )
    
def resample(x, xp, yp, xedges=False, xpedges=False):
    """
    IN PROGRESS.  Resample a spectrum to a new binning using PixelSpline

    1 <= x.ndim <= xp.ndim <= yp.ndim <= 2
    """

    assert 1 <= x.ndim
    assert x.ndim <= xp.ndim
    assert xp.ndim <= yp.ndim
    assert yp.ndim <= 2

    input_edges = xp if xpedges else pixspline.cen2bound(xp)
    ys = pixspline.PixelSpline(input_edges, yp)
    
    edges = x if xedges else pixspline.cen2bound(x)
    return ys.resample(edges)
    
        
    
    
    
    
    
    
    
    
    
    
