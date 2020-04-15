"""
specter.util.traceset
=====================

Handle sets of Legendre coefficients
"""

from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import os
import numpy as np
from numpy.polynomial.legendre import legfit
import numbers
from specter.util import legval_numba

class TraceSet(object):
    def __init__(self, coeff, domain=[-1,1]):
        """
        coeff[nspec, deg+1]
        map dependent variable within domain to [-1,1]
        """
        self._coeff = coeff.copy()
        self._xmin = domain[0]
        self._xmax = domain[1]

    @property
    def ntrace(self):
        return self._coeff.shape[0]

    def _xnorm(self, x):
        #return 2.0 * (x - self._xmin) / (self._xmax - self._xmin) - 1.0
        #implement oleksandr-pavlyk's fix from the memory branch
        return (x - self._xmin) * (2.0 / (self._xmax - self._xmin)) - 1.0

    def eval(self, ispec, x):
        # TODO: this function could be a little more robust/elegant
        xx = np.array(self._xnorm(x))
        #wavelength input may be a scalar
        scalar_input = np.isscalar(x)
        #ispec input may be in the form (min, max)
        tuple_input = False
        if type(ispec) is tuple:
            tuple_input = True
            spec_min, spec_max = ispec
            #in this case ispec is coming from the cached branch of xypix, which means
            #that we need to use nspec instead of ispec due to re-indexing
            nspec = spec_max - spec_min
            nwave = len(x)

        #numba requires f8 or smaller
        if tuple_input is False:
            cc_numba = self._coeff[ispec].astype(np.float64, copy=False)
        else:
            cc_numba = self._coeff[spec_min:spec_max].astype(np.float64, copy=False)

        #use numba version of legval if possible
        if isinstance(ispec, numbers.Integral):
            results = legval_numba(xx, cc_numba)
            if scalar_input:
                return results[0]
            else:
                return results
        #handle tuple case from _xypix cached branch
        elif tuple_input:
            #compute and store the values
            y=np.zeros([nspec, nwave])
            for i in range(nspec):
                y[i,:]=legval_numba(xx, cc_numba[i])
            return y

        else:
            #ispec may sometimes be None
            if ispec is None:
                ispec = list(range(self._coeff.shape[0]))

            #use numba version of legval if possible
            y=[]
            for i in ispec:
                cc_i = self._coeff[i].astype(np.float64, copy=False)
                if scalar_input:
                    y.append(legval_numba(xx, cc_i)[0])
                else:
                    y.append(legval_numba(xx, cc_i))

            return np.array(y)


    # def __call__(self, ispec, x):
    #     return self.eval(ispec, x)

    def invert(self, domain=None, deg=None):
        """
        Return a traceset modeling x vs. y instead of y vs. x
        """
        ytmp = self.eval(None, np.array([self._xmin, self._xmax]))
        ymin = np.min(ytmp)
        ymax = np.max(ytmp)
        x = np.linspace(self._xmin, self._xmax, 1000)
        if deg is None:
            deg = self._coeff.shape[1]+2

        c = np.zeros((self.ntrace, deg+1))
        for i in range(self.ntrace):
            y = self.eval(i, x)
            #yy = 2.0 * (y-ymin) / (ymax-ymin) - 1.0
            #implement oleksandr-pavlyk's fix from the memory branch
            yy = (y-ymin) * (2.0 / (ymax-ymin)) - 1.0
            c[i] = legfit(yy, x, deg)

        return TraceSet(c, domain=(ymin, ymax))

def fit_traces(x, yy, deg=5, domain=None):
    """
    returns TraceSet object modeling y[i] vs. x

    Args:
        x : 1D array
        y : 2D array[nspec, nx]
        deg : optional Legendre degree
    """
    nspec, nx = yy.shape
    assert len(x) == nx, "y.shape[1] ({}) != len(x) ({})".format(nx, len(x))
    assert np.all(np.diff(x) > 0), "x not monotonically increasing"

    if domain is None:
        xmin, xmax = x[0], x[-1]
    else:
        xmin, xmax = domain

    c = np.zeros((nspec, deg+1))
    #implement oleksandr-pavlyk's fix from the memory branch
    xx = x - xmin
    xx *= 2.0/(xmax-xmin)
    xx -= 1.0
    for i in range(nspec):
        #xx = 2.0 * (x-xmin) / (xmax-xmin) - 1.0
        c[i] = legfit(xx, yy[i], deg)

    return TraceSet(c, [xmin, xmax])
