"""
Handle sets of Legendre coefficients
"""

import sys
import os
import numpy as N
from numpy.polynomial.legendre import legfit, legval

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
        if not isinstance(x, (int,float,N.ndarray)):
            x = N.array(x)
        return 2.0 * (x - self._xmin) / (self._xmax - self._xmin) - 1.0
        
    def eval(self, ispec, x):
        xx = self._xnorm(x)
        
        if isinstance(ispec, int):
            return legval(xx, self._coeff[ispec])
        else:
            if ispec is None:
                ispec = range(self._coeff.shape[0])
            
            y = [legval(xx, self._coeff[i]) for i in ispec]
            return N.array(y)
            
    # def __call__(self, ispec, x):
    #     return self.eval(ispec, x)
            
    def invert(self, domain=None, deg=None):
        """
        Return a traceset modeling x vs. y instead of y vs. x
        """
        ytmp = self.eval(None, (self._xmin, self._xmax))
        ymin = N.min(ytmp)
        ymax = N.max(ytmp)
        x = N.linspace(self._xmin, self._xmax, 1000)
        if deg is None:
            deg = self._coeff.shape[1]+2
            
        c = N.zeros((self.ntrace, deg+1))
        for i in range(self.ntrace):
            y = self.eval(i, x)
            yy = 2.0 * (y-ymin) / (ymax-ymin) - 1.0
            c[i] = legfit(yy, x, deg)
            
        return TraceSet(c, domain=(ymin, ymax))
            
        
        
    
