"""
Handle sets of Legendre coefficients
"""

from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import os
import numpy as np
from numpy.polynomial.legendre import legfit, legval
import numbers

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
        
    def xnorm(self, x):
        '''
        Convert x -> [-1,1]
        '''
        if not isinstance(x, (numbers.Real, np.ndarray)):
            x = np.array(x)
        return (x - self._xmin) * (2.0 / (self._xmax - self._xmin)) - 1.0
        
    def eval(self, ispec, x=None, xnorm=None):
        '''
        Evaluate traceset at values x
        if x1 is provided, it is the "reduced x" converted to domain [-1,1]
        '''
        if xnorm is None:
            xnorm = self.xnorm(x)
        
        if isinstance(ispec, numbers.Integral):
            return legval(xnorm, self._coeff[ispec])
        else:
            if ispec is None:
                ispec = list(range(self._coeff.shape[0]))

            #y = [legval(xnorm, self._coeff[i]) for i in ispec]
            #return np.array(y)
            return legval(xnorm, self._coeff[ispec].T)
            
    # def __call__(self, ispec, x):
    #     return self.eval(ispec, x)
            
    def invert(self, domain=None, deg=None):
        """
        Return a traceset modeling x vs. y instead of y vs. x
        """
        ytmp = self.eval(None, (self._xmin, self._xmax))
        ymin = np.min(ytmp)
        ymax = np.max(ytmp)
        x = np.linspace(self._xmin, self._xmax, 1000)
        if deg is None:
            deg = self._coeff.shape[1]+2
            
        c = np.zeros((self.ntrace, deg+1))
        for i in range(self.ntrace):
            y = self.eval(i, x)
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
    # xx = 2.0 * (x-xmin) / (xmax-xmin) - 1.0
    xx = x - xmin
    xx *= 2.0/(xmax-xmin)
    xx -= 1.0
    for i in range(nspec):
        c[i] = legfit(xx, yy[i], deg)

    return TraceSet(c, [xmin, xmax])
    
    
    
        
    
