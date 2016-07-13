#!/usr/bin/env python

"""
Unit tests for PSF classes.
"""

from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import os
import numpy as np
from numpy.polynomial import legendre
from specter import util
import unittest
from pkg_resources import resource_filename
from specter.psf import load_psf

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
        x = np.linspace(0, 2, 20)
        y = np.sin(x)
        #- Check integration over full range
        edges = (x[0], x[-1])
        self.assertTrue(np.trapz(y, x) == util.trapz(edges, x, y)[0])

        #- Check that integrations add up to each other
        lo = 0.5*(x[0] + x[1])
        hi = 0.5*(x[-2] + x[-1])
        mid = 0.5*(lo + hi)
        edges = np.linspace(lo, hi, 30)
        blat = util.trapz(edges, x, y)
        foo = util.trapz([lo, hi], x, y)
        bar = util.trapz([lo, mid, hi], x, y)
        self.assertAlmostEqual(np.sum(blat), foo[0])
        self.assertAlmostEqual(np.sum(blat), np.sum(bar))

        #- Edges outside of range shouldn't crash
        blat = util.trapz([-1,0,1,2], x, y)

    def test_sincshift(self):
        a = np.zeros((3,5))
        self.assertEqual(a.shape, util.sincshift(a, 0.1, 0.0).shape)
        self.assertEqual(a.shape, util.sincshift(a, 0.0, 0.1).shape)
        self.assertEqual(a.shape, util.sincshift(a, 0.1, 0.1).shape)
        self.assertEqual(a.shape, util.sincshift2d(a, 0.1, 0.0).shape)
        self.assertEqual(a.shape, util.sincshift2d(a, 0.0, 0.1).shape)
        self.assertEqual(a.shape, util.sincshift2d(a, 0.1, 0.1).shape)

    def test_psfbias(self):
        psf = load_psf(resource_filename('specter.test', "t/psf-pix.fits"))
        wmid = 0.5*(psf.wmin+psf.wmax)
        ww = np.linspace(wmid-10, wmid+10)
        phot = np.ones(len(ww))
        phot[10] = 20
        bias = util.psfbias(psf, psf, ww, phot)
        absbias, R = util.psfabsbias(psf, psf, ww, phot)

    def test_resample(self):
        '''test resample; coverage only (not actual functionality)'''
        x = np.linspace(0, 2*np.pi, 20)
        y = np.sin(x) + 1
        newx = np.linspace(1, 6, 10)
        newy = util.resample(newx, x, y)

    def test_timeit(self):
        import time
        t0 = util.util._timeit()
        time.sleep(0.1)
        dt = util.util._timeit()
        self.assertGreater(dt, 0.1)
        self.assertLess(dt, 0.11)

    # def test_rebin(self):
    #     x = np.arange(25)
    #     y = np.random.uniform(0.0, 5.0, size=len(x))
    #     y[0] = y[-1] = 0.0  #- avoid edge bin effects for comparison
    #     xx = np.arange(0, len(x), 1.3)
    #     yy = util.rebin(xx, x, y)
    #
    #     s1 = np.sum(y * np.gradient(x))
    #     s2 = np.sum(yy * np.gradient(xx))
    #     self.assertAlmostEqual(s1, s2)
    #
    #     edges = np.linspace(0, len(x))
    #     zz = util.rebin(edges, x, y, edges=True)
    #     s3 = np.sum(yy * np.gradient(xx))
    #     self.assertEqual(len(zz), len(edges)-1)
    #     self.assertAlmostEqual(s1, s3)
    #
    #     with self.assertRaises(ValueError):
    #         zz = util.rebin(xx, x[-1::-1], y)

if __name__ == '__main__':
    unittest.main()
