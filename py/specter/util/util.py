"""
Utility functions and classes for specter

Stephen Bailey
Fall 2012
"""

from __future__ import absolute_import, division, print_function, unicode_literals

import os
import math
import numpy as np
import scipy.signal
from scipy.special import legendre
from scipy.sparse import spdiags
from scipy.signal import convolve, convolve2d
from specter.util import pixspline
from time import time

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
        self.x = np.array(x)
        self.y = np.array(y)
        self.data = np.array(data)

    def __call__(self, x, y):
        """
        Evaluate data at (x,y)
        """
        #- TODO: compare speed to solution at
        #- http://stackoverflow.com/questions/12729228/simple-efficient-bilinear-interpolation-of-images-in-numpy-and-python

        #- Find where we are in grid
        #- clip to 1 because we will use i and i-1
        #- clip to len(x)-1 to allow extrapolation beyond grid boundary
        ix = np.searchsorted(self.x, x).clip(1, len(self.x)-1)
        iy = np.searchsorted(self.y, y).clip(1, len(self.y)-1)

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
        xx = (x+dx)*np.pi  #- cache shifted array for 30% faster evals
        return np.exp( -(xx/(dampfac*np.pi))**2 ) * np.sin(xx) / xx
    else:
        xx = np.zeros(len(x))
        xx[len(x)//2] = 1.0
        return xx

#- Implementation note: the typical PSF image is 15x15.
#- fftconvolve is not faster than convolve for images this small
def sincshift(image, dx, dy, sincrad=10, dampfac=3.25):
    """
    Return image shifted by dx, dy using sinc interpolation.

    For speed, do each dimension independently which can introduce edge
    effects.  Also see sincshift2d().
    """
    s = np.arange(-sincrad, sincrad+1.0)
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
    s = np.arange(-sincrad, sincrad+1.0)
    sincx = _sincfunc(s, -dx, dampfac=dampfac)
    sincy = _sincfunc(s, -dy, dampfac=dampfac)
    kernel = np.outer(sincy, sincx)
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
    edges = np.concatenate((x-0.5, x[-1:]+0.5))
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
    x = np.linalg.lstsq(iCov, y)[0]
    return x, iCov

def trapz(edges, xp, yp):
    """
    Perform trapezoidal integration between edges using sampled function
    yp vs. xp.  Returns array of length len(edges)-1.

    Input xp array must be sorted in ascending order.

    See also numpy.trapz, which integrates a single array
    """
    if np.any(np.diff(xp) < 0.0):
        raise ValueError("Input x must be sorted in increasing order")

    if len(xp) != len(yp):
        raise ValueError("xp and yp must have same length")

    yedge = np.interp(edges, xp, yp)
    result = np.zeros(len(edges)-1)
    iedge = np.searchsorted(xp, edges)
    for i in range(len(edges)-1):
        ilo, ihi = iedge[i], iedge[i+1]
        xx = np.concatenate( (edges[i:i+1], xp[ilo:ihi], edges[i+1:i+2]) )
        yy = np.concatenate( (yedge[i:i+1], yp[ilo:ihi], yedge[i+1:i+2]) )
        result[i] = np.trapz(yy, xx)

    return result

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

#- Faster versions than np.outer, which has to do type checking and raveling
try:
    #- ~3x faster if numba is installed
    if 'NUMBA_DISABLE_JIT' in os.environ:
        raise ImportError
    import numba
    @numba.jit
    def outer(x, y, out):
        for i in range(len(x)):
            for j in range(len(y)):
                out[i,j] = x[i] * y[j]
        return out

except ImportError:
    #- 1.5x faster otherwise
    def outer(x, y, out):
        return np.multiply(x[:, None], y[None, :], out)


# Much faster than numpy.polynomial.legendre.legval, but doesn't work with scalars
import numba
@numba.jit(nopython=True,cache=False)
def legval_numba(x, c):
    nd=len(c)
    ndd=nd
    xlen = x.size
    c0=c[-2]*np.ones(xlen)
    c1=c[-1]*np.ones(xlen)
    for i in range(3, ndd + 1):
        tmp = c0
        nd = nd - 1
        nd_inv = 1/nd
        c0 = c[-i] - (c1*(nd - 1))*nd_inv
        c1 = tmp + (c1*x*(2*nd - 1))*nd_inv
    return c0 + c1*x


@numba.jit(nopython=True, cache=False)
def custom_hermitenorm(n, u):
    """
    Custom implementation of scipy.special.hermitenorm to enable jit-compiling
        with Numba (which as of 10/2018 does not support scipy). This functionality
        is equivalent to:
        fn = scipy.special.hermitenorm(n)
        return fn(u)
        with the exception that scalar values of u are not supported.
    Inputs:
        n: the degree of the hermite polynomial
        u: (requires array) points at which the polynomial will be evaulated.
    Outputs:
        res: the value of the hermite polynomial at array points(u)
    """

    #below is (mostly) cut and paste from scipy orthogonal_eval.pxd
    #some modifications have been made to operate on an array
    #rather than a single value (as in the original version)
    res=np.zeros(len(u))
    if n < 0:
        return (0.0)*np.ones(len(u))
    elif n == 0:
        return (1.0)*np.ones(len(u))
    elif n == 1:
        return u
    else:
        y3 = 0.0
        y2 = 1.0
        for i,x in enumerate(u):
            for k in range(n, 1, -1):
                y1 = x*y2-k*y3
                y3 = y2
                y2 = y1
            res[i]=x*y2-y3
            #have to reset before the next iteration
            y3 = 0.0
            y2 = 1.0
        return res


@numba.jit(nopython=True, cache=False)
def custom_erf(y):
    """Custom implementation of :func:`scipy.special.erf` to enable jit-compiling
    with Numba (which as of 10/2018 does not support scipy). This functionality is equilvalent to::

        scipy.special.erf(y)

    with the exception that scalar values of y are not supported.

    Parameters
    ----------
    y : array-like
        Points at which the error function will be evaluated.

    Returns
    -------
    array-like
        The value of the error function at points in array `y`.

    Notes
    -----
    This function has been translated from the original fortran function
    to Python. The original scipy erf function can be found at:
    https://github.com/scipy/scipy/blob/8dba340293fe20e62e173bdf2c10ae208286692f/scipy/special/cdflib/erf.f
    Note that this new function introduces a small amount of machine-precision numerical error
    as compared to the original scipy function.
    """

    #have to define a ton of constants
    c=0.564189583547756E0
    ###
    a1=0.771058495001320E-04
    a2=-0.133733772997339E-02
    a3=0.323076579225834E-01
    a4=0.479137145607681E-01
    a5=0.128379167095513E+00
    ###
    b1=0.301048631703895E-02
    b2=0.538971687740286E-01
    b3=0.375795757275549E+00
    ###
    p1=-1.36864857382717E-07
    p2=5.64195517478974E-01
    p3=7.21175825088309E+00
    p4=4.31622272220567E+01
    p5=1.52989285046940E+02
    p6=3.39320816734344E+02
    p7=4.51918953711873E+02
    p8=3.00459261020162E+02
    ###
    q1=1.00000000000000E+00
    q2=1.27827273196294E+01
    q3=7.70001529352295E+01
    q4=2.77585444743988E+02
    q5=6.38980264465631E+02
    q6=9.31354094850610E+02
    q7=7.90950925327898E+02
    q8=3.00459260956983E+02
    ###
    r1=2.10144126479064E+00
    r2=2.62370141675169E+01
    r3=2.13688200555087E+01
    r4=4.65807828718470E+00
    r5=2.82094791773523E-01
    ###
    s1=9.41537750555460E+01
    s2=1.87114811799590E+02
    s3=9.90191814623914E+01
    s4=1.80124575948747E+01
    ###
    #end of constants

    #the orig version is meant for a single point
    #need to modify to work on an array
    erf = np.zeros(len(y))
    for i,x in enumerate(y):
        ax=abs(x)
        #change gotos into something sensible
        if ax <= 0.5E0:
            t=x*x
            top = ((((a1*t+a2)*t+a3)*t+a4)*t+a5) + 1.0E0
            bot = ((b1*t+b2)*t+b3)*t + 1.0E0
            erf[i] = x * (top/bot)
        elif 0.5E0 < ax <= 4.0E0:
            top = ((((((p1*ax+p2)*ax+p3)*ax+p4)*ax+p5)*ax+p6)*ax + p7)*ax + p8
            bot = ((((((q1*ax+q2)*ax+q3)*ax+q4)*ax+q5)*ax+q6)*ax + q7)*ax + q8
            val = 0.5E0 + (0.5E0 - np.exp(-x*x)*top/bot)
            if x < 0.0E0:
                erf[i] = -val
            else:
                erf[i] = val
        elif 4.0E0 < ax < 5.8E0:
            x2 = x*x
            t = 1.0E0/x2
            top = (((r1*t+r2)*t+r3)*t+r4)*t + r5
            bot = (((s1*t+s2)*t+s3)*t+s4)*t + 1.0E0
            val = (c-top/ (x2*bot)) / ax
            val = 0.5E0 + (0.5E0 - np.exp(-x2)*val)
            if x < 0.0E0:
                erf[i] = -val
            else:
                erf[i] = val
        elif ax >= 5.8E0:
            #choose the sign
            if x < 0.0E0:
                erf[i] = -1.0E0
            else:
                erf[i] = 1.0E0

    return erf
