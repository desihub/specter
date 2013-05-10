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

    def test_trapz(self):
        x = N.linspace(0, 2, 20)
        y = N.sin(x)
        #- Check integration over full range
        edges = (x[0], x[-1])
        self.assertTrue(N.trapz(y, x) == util.trapz(edges, x, y)[0])
        
        #- Check that integrations add up to each other
        lo = 0.5*(x[0] + x[1])
        hi = 0.5*(x[-2] + x[-1])
        mid = 0.5*(lo + hi)
        edges = N.linspace(lo, hi, 30)
        blat = util.trapz(edges, x, y)
        foo = util.trapz([lo, hi], x, y)
        bar = util.trapz([lo, mid, hi], x, y)
        self.assertTrue(N.sum(blat) == foo[0])
        self.assertTrue(N.sum(blat) == N.sum(bar))
        
        #- Edges outside of range shouldn't crash
        blat = util.trapz([-1,0,1,2], x, y)
        
    # def test_rebin(self):
    #     x = N.arange(25)
    #     y = N.random.uniform(0.0, 5.0, size=len(x))
    #     y[0] = y[-1] = 0.0  #- avoid edge bin effects for comparison
    #     xx = N.arange(0, len(x), 1.3)
    #     yy = util.rebin(xx, x, y)
    #     
    #     s1 = N.sum(y * N.gradient(x))
    #     s2 = N.sum(yy * N.gradient(xx))
    #     self.assertAlmostEqual(s1, s2)
    #     
    #     edges = N.linspace(0, len(x))
    #     zz = util.rebin(edges, x, y, edges=True)
    #     s3 = N.sum(yy * N.gradient(xx))
    #     self.assertEqual(len(zz), len(edges)-1)
    #     self.assertAlmostEqual(s1, s3)
    #     
    #     with self.assertRaises(ValueError):
    #         zz = util.rebin(xx, x[-1::-1], y)

if __name__ == '__main__':
    unittest.main()           
