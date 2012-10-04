"""
Utility functions and classes for specter

Stephen Bailey
Fall 2012
"""

import numpy as N
import scipy.signal

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
        #- Find where we are in grid
        #- clip to 1 because we will use i and i-1
        #- clip to len(x)-1 to allow extrapolation beyond grid boundary
        ix = N.searchsorted(self.x, x).clip(1, len(self.x)-1)
        iy = N.searchsorted(self.y, y).clip(1, len(self.y)-1)
        
        #- Interpolation distances from points
        dx = (x - self.x[ix-1]) / (self.x[ix] - self.x[ix-1])
        dy = (y - self.y[iy-1]) / (self.y[iy] - self.y[iy-1])
        
        #- Interpolate, allowing x and/or y to be multi-dimensional
        data1 = (self.data[ix-1,iy-1].T*(1-dx) + self.data[ix,iy-1].T*dx).T
        data2 = (self.data[ix-1,iy].T*(1-dx) + self.data[ix,iy].T*dx).T
        dataxy = (data1.T*(1-dy) + data2.T*dy).T
        
        return dataxy
        
def psfbias(pix1, pix2):
    """
    Return bias from extracting true PSF with pix1 using PSF with pix2
    """
    return 1.0 - N.sum(pix1*pix2) / N.sum(pix1*pix1)

def rebin(pix, n):
    """
    rebin 2D array pix into bins of size n x n
    
    New binsize must be evenly divisible into original pix image
    """
    assert pix.shape[0] % n == 0
    assert pix.shape[1] % n == 0
    
    s = pix.shape[0]//n, n, pix.shape[1]//n, n
    return pix.reshape(s).sum(-1).sum(1)

    
#- Utility functions for sinc shifting pixelated PSFs
def _sincfunc(x, dx, dampfac=3.25):
    """sinc helper function for sincshift()"""
    if dx != 0.0:
        return N.exp( -((x+dx)/dampfac)**2 ) * N.sin( N.pi*(x+dx) ) / (N.pi * (x+dx))
    else:
        xx = N.zeros(len(x))
        xx[len(x)/2] = 1.0
        return xx

def sincshift(image, dx, dy, sincrad=10, dampfac=3.25):
    """
    Return image shifted by dx, dy using sinc interpolation
    """
    s = N.arange(-sincrad, sincrad+1)
    sincx = _sincfunc(s, -dx, dampfac=dampfac)

    #- If we're shifting just in x, do faster 1D convolution with wraparound
    #- WARNING: can introduce edge effects if image isn't nearly 0 at edge
    if abs(dy) < 1e-6:
        newimage = scipy.signal.convolve(image.ravel(), sincx, mode='same')
        return newimage.reshape(image.shape)

    sincy = _sincfunc(s, -dy, dampfac=dampfac)
    kernel = N.outer(sincy, sincx)
    newimage = scipy.signal.convolve2d(image, kernel, mode='same')
    return newimage

#- Integrals
# def trapz(x, xp, yp):
#     """
#     Given a function yp[] sampled at values xp[], return an array of
#     size len(xedges)-1 integating yp(xp) between the values given in x
#     using trapezoidal integration.
#     """
#     raise NotImplementedError
