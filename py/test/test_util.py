#!/usr/bin/env python

"""
Unit tests for PSF classes.
"""

import sys
import os
import numpy as N
from numpy.polynomial import legendre
from specter import util
import unittest

class TestUtil(unittest.TestCase):
    """
    Test functions within specter.util
    """
    
    def test_gaussX(self):
        #- Gaussian integration
        self.assertTrue(util.gaussint(-100) == 0.0)
        self.assertTrue(util.gaussint(0.0) == 0.5)
        self.assertTrue(util.gaussint(+100) == 1.0)
        
        for x in (-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 10.0):
            self.assertTrue(util.gaussint(x, mean=x, sigma=1.0) == 0.5)
            self.assertTrue(util.gaussint(x, mean=x, sigma=2.0) == 0.5)
            if x>0:
                self.assertTrue(util.gaussint(x, sigma=x) == 0.84134474606854293)
                self.assertTrue(util.gaussint(-x, sigma=x) == 0.15865525393145707)

    def test_resample(self):
        x = N.arange(25)
        y = N.random.uniform(0.0, 5.0, size=len(x))
        y[0] = y[-1] = 0.0  #- avoid edge bin effects for comparison
        xx = N.arange(0, len(x), 1.3)
        yy = util.resample(x, y, xnew=xx)
        
        s1 = N.sum(y * N.gradient(x))
        s2 = N.sum(yy * N.gradient(xx))
        self.assertAlmostEqual(s1, s2)
        
        edges = N.linspace(0, len(x))
        zz = util.resample(x, y, edges=edges)
        s3 = N.sum(yy * N.gradient(xx))
        self.assertEqual(len(zz), len(edges)-1)
        self.assertAlmostEqual(s1, s3)
        
        with self.assertRaises(ValueError):
            zz = util.resample(x[-1::-1], y, xnew=xx)        

if __name__ == '__main__':
    unittest.main()           
