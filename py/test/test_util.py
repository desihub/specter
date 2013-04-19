#!/usr/bin/env python

"""
Unit tests for PSF classes.
"""

import sys
import os
import numpy as N
from numpy.polynomial import legendre
import specter.util
import unittest

class TestUtil(unittest.TestCase):
    """
    Test functions within specter.util
    """
    
    def test_LegendreFit(self):
        x = N.arange(10.0)
        xx = 2 * (x-x.min()) / (x.max() - x.min()) - 1.0
        c = N.array( (1.0, 1.5, -0.1) )
        y = legendre.legval(xx, c)
        
        #- Test input points == output fit at same points
        fit = specter.util.LegendreFit(x, y)
        self.assertTrue(N.all(N.round(y-fit(x), 12) == 0.0))
        
        #- Test coefficients
        nc = len(c)
        self.assertTrue(N.all(N.round(c - fit.coeff[0:nc], 12) == 0.0))
        self.assertTrue(N.all(N.abs(fit.coeff[nc:]) < 1e-12))
        
        #- Test fit for other ranges
        x = N.linspace(0.0, 10.0)
        y = fit(x)
        self.assertTrue(len(x) == len(y))  #- not very stringent...
        
    def test_gaussX(self):
        #- Gaussian integration
        self.assertTrue(specter.util.gaussint(-100) == 0.0)
        self.assertTrue(specter.util.gaussint(0.0) == 0.5)
        self.assertTrue(specter.util.gaussint(+100) == 1.0)
        
        for x in (-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 10.0):
            self.assertTrue(specter.util.gaussint(x, mean=x, sigma=1.0) == 0.5)
            self.assertTrue(specter.util.gaussint(x, mean=x, sigma=2.0) == 0.5)
            if x>0:
                self.assertTrue(specter.util.gaussint(x, sigma=x) == 0.84134474606854293)
                self.assertTrue(specter.util.gaussint(-x, sigma=x) == 0.15865525393145707)
            

if __name__ == '__main__':
    unittest.main()           
