#!/usr/bin/env python

"""
Test PSF classes
"""

import sys
import os
import numpy as N
from specter.psf import load_psf

import unittest

class TestPSF(unittest.TestCase):
    # def setUp(self):
    #     self.indir = os.environ['SPECTER_DIR'] + '/test/data'
    #     self.psf = load_psf(self.indir + "/psf-pix.fits")

    def wrap_wave_test(self, fn):
        #- Unspecified y -> array with nwave elements
        tmp = fn(0, None)
        self.assertTrue(len(tmp) == self.psf.nwave)
        
        #- Specified fiber and y -> float
        tmp = fn(0, y=0)
        self.assertTrue(isinstance(tmp, float))

        #- Arbitrary
        yy = N.linspace(0, self.psf.nwave)
        tmp = fn(0, y=yy)
        self.assertTrue(len(tmp) == len(yy))
        
        #- Full loglam grid
        self.assertTrue(fn().shape == (self.psf.nspec, self.psf.nwave))
        
        #- ispec >= nspec should raise an error
        with self.assertRaises(IndexError):
            fn(self.psf.nspec)

    def test_loglam(self):
        self.wrap_wave_test(self.psf.loglam)

    def test_wavelength(self):
        self.wrap_wave_test(self.psf.wavelength)
        
    def test_dimensions(self):
        self.assertTrue(self.psf.npix_x > 0)
        self.assertTrue(self.psf.npix_y > 0)
        self.assertTrue(self.psf.nspec > 0)
        self.assertTrue(self.psf.nwave > 0)
        
    ### @unittest.skip("BROKEN")
    def test_pix(self):
        ww = self.psf.wavelength()
        wtest = list()
        wtest.append(N.min(ww[:, 0]))
        wtest.append(N.max(ww[:, 0]))
        wtest.append(N.min(ww[:, -1]))
        wtest.append(N.max(ww[:, -1]))
        wtest.append(N.mean(wtest))
        
        for i in (0, self.psf.nspec/2, self.psf.nspec-1):
            for w in wtest:
                pix = self.psf.pix(i, w)
                self.assertEquals(pix.ndim, 2)  

    ### @unittest.skip("BROKEN")
    def test_xypix(self):
        ww = self.psf.wavelength()
        wtest = list()
        wtest.append(N.min(ww[:, 0]))
        wtest.append(N.max(ww[:, 0]))
        wtest.append(N.min(ww[:, -1]))
        wtest.append(N.max(ww[:, -1]))
        wtest.append(N.mean(wtest))
        
        for i in (0, self.psf.nspec/2, self.psf.nspec-1):
            for w in wtest:
                xx, yy, pix = self.psf.xypix(self.psf.nspec/2, w)
                shape = (yy.stop-yy.start, xx.stop-xx.start)
                msg = "%s != %s at (i=%d, w=%.1f)" % (str(pix.shape), str(shape), i, w)
                self.assertEqual(pix.shape, shape, msg)
                
    def test_xypix_range(self):
        w = N.mean(self.psf.wavelength())
        i = self.psf.nspec/2
        x0, y0, pix = self.psf.xypix(i, w)
        xx, yy, pix = self.psf.xypix(i, w, xmin=x0.start)
        self.assertTrue(xx.start == 0)
        self.assertTrue(yy.start == y0.start)
        xx, yy, pix = self.psf.xypix(i, w, ymin=y0.start)
        self.assertTrue(xx.start == x0.start)
        self.assertTrue(yy.start == 0)
        xx, yy, pix = self.psf.xypix(i, w, xmax=x0.stop-1)
        self.assertTrue(xx.start == x0.start)
        self.assertTrue(xx.stop == x0.stop-1)
        xx, yy, pix = self.psf.xypix(i, w, ymax=y0.stop-1)
        self.assertTrue(yy.start == y0.start)
        self.assertTrue(yy.stop == y0.stop-1)
        xx, yy, pix = self.psf.xypix(i, w, xmin=x0.start, ymin=y0.start)
        self.assertTrue(xx.start == 0)
        self.assertTrue(yy.start == 0)
        
        
    def test_project11(self):
        ww = self.psf.wavelength(0)[0:10]
        phot = N.random.uniform(0,100,len(ww))
        img = self.psf.project(phot, ww, verbose=False)
        self.assertEquals(img.shape, (self.psf.npix_y, self.psf.npix_x))
    
    def test_project12(self):
        ww = self.psf.wavelength(0)[0:10]
        phot = N.random.uniform(0,100,len(ww))
        phot = N.tile(phot, 5).reshape(5, len(ww))
        img = self.psf.project(phot, ww, verbose=False)
        self.assertEquals(img.shape, (self.psf.npix_y, self.psf.npix_x))

    ### @unittest.skip("BROKEN")
    def test_project22(self):
        ww = self.psf.wavelength(0)[0:10]
        ww = N.tile(ww, 5).reshape(5, len(ww))
        phot = N.random.uniform(0,100,len(ww))
        phot = N.tile(phot, 5).reshape(5, len(ww))
        img = self.psf.project(phot, ww, verbose=False)
        self.assertEquals(img.shape, (self.psf.npix_y, self.psf.npix_x))
    
    def test_project_specmin(self):
        ww = self.psf.wavelength(0)[0:10]
        phot = N.random.uniform(0,100,len(ww))
        img = self.psf.project(phot, ww, specmin=1, verbose=False)
        self.assertEquals(img.shape, (self.psf.npix_y, self.psf.npix_x))
        with self.assertRaises(IndexError):
            i = self.psf.nspec
            img = self.psf.project(phot, ww, specmin=i, verbose=False)

    #- Needs psf.xyrange to work first
    @unittest.skip("Skipping")
    def test_project_xyrange(self):
        self.assertTrue(True)

    def test_shift_xy(self):
        dx, dy = 0.1, 0.2
        x0 = self.psf.x(0).copy()
        y0 = self.psf.y(0).copy()
        self.psf.shift_xy(dx, dy)
        self.assertTrue(N.all(self.psf.x(0) == x0+dx))
        self.assertTrue(N.all(self.psf.y(0) == y0+dy))
    
    def wrap_xy_test(self, xyfn):
        w = self.psf.wavelength(0)[0:10]
        x = xyfn(0, w)
        self.assertEqual(x.shape, w.shape)
    
    def test_x(self):
        #- Grid of x positions
        x = self.psf.x()
        self.assertEqual(x.ndim, 2)
        self.assertEqual(x.shape[0], self.psf.nspec)
        self.assertEqual(x.shape[1], self.psf.nwave)

        #- x for fiber 0
        x = self.psf.x(0)
        self.assertEqual(x.ndim, 1)
        self.assertEqual(len(x), self.psf.nwave)
        
        #- x for fiber 0 at a specific wavelength
        w = N.mean(self.psf.wavelength(0))
        x = self.psf.x(0, w)
        self.assertTrue(isinstance(x, float))
        
        #- x for all fibers at a fixed wavelength
        x = self.psf.x(None, w)
        self.assertEqual(x.ndim, 1)
        self.assertEqual(len(x), self.psf.nspec)
        
        #- x for fiber 0 at a range of wavelengths
        w = self.psf.wavelength(0)[0:10]
        x = self.psf.x(0, w)
        self.assertEqual(x.shape, w.shape)

        #- x for all fibers at a range of wavelengths
        x = self.psf.x(None, w)
        self.assertEqual(x.shape, (self.psf.nspec, len(w)))

    ### @unittest.skip("BROKEN")        
    def test_y(self):
        #- Grid of y positions
        y = self.psf.y()
        self.assertEqual(y.ndim, 2)
        self.assertEqual(y.shape[0], self.psf.nspec)
        self.assertEqual(y.shape[1], self.psf.nwave)

        #- y for fiber 0
        y = self.psf.y(0)
        self.assertEqual(y.ndim, 1)
        self.assertEqual(len(y), self.psf.nwave)
        
        #- y for fiber 0 at a specific wavelength
        w = N.mean(self.psf.wavelength(0))
        y = self.psf.y(0, w)
        #- FAILS: returns 0dim array !?!
        self.assertTrue(isinstance(y, float))
        
        #- y for all fibers at a fixed wavelength
        y = self.psf.y(None, w)
        self.assertEqual(y.ndim, 1)
        self.assertEqual(len(y), self.psf.nspec)
        
        #- y for fiber 0 at a range of wavelengths
        w = self.psf.wavelength(0)[0:10]
        y = self.psf.y(0, w)
        self.assertEqual(y.shape, w.shape)

        #- y for all fibers at a range of wavelengths
        y = self.psf.y(None, w)
        self.assertEqual(y.shape, (self.psf.nspec, len(w)))
        
    def test_xy(self):
        x = self.psf.x(0)
        y = self.psf.y(0)
        xy = self.psf.xy(0)
        self.assertTrue(N.all(xy[0] == x))
        self.assertTrue(N.all(xy[1] == y))
        
    def test_xyw(self):
        x = self.psf.x(0)
        y = self.psf.y(0)
        w = self.psf.wavelength(0)
        xyw = self.psf.xyw(0)
        self.assertTrue(N.all(xyw[0] == x))
        self.assertTrue(N.all(xyw[1] == y))
        self.assertTrue(N.all(xyw[2] == w))
    
    #- BROKEN
    ### @unittest.skip("BROKEN")
    def test_xyrange(self):
        ww = self.psf.wavelength(0)[0:10]
        wave_range = [N.min(ww), N.max(ww)]
        spec_range = [0,10]
        xmin,xmax,ymin,ymax = self.psf.xyrange(spec_range, wave_range)
        # print xmin,xmax,ymin,ymax
        for w in wave_range:
            for i in spec_range:
                xx, yy, pix = self.psf.xypix(i, w)
                self.assertGreaterEqual(xx.start, xmin)
                self.assertLessEqual(xx.stop, xmax)
                self.assertGreaterEqual(yy.start, ymin)
                self.assertLessEqual(yy.stop, ymax)
            
class TestPixPSF(TestPSF):
    def setUp(self):
        self.indir = os.environ['SPECTER_DIR'] + '/test/data'
        self.psf = load_psf(self.indir + "/psf-pix.fits")

class TestSpotPSF(TestPSF):
    def setUp(self):
        self.indir = os.environ['SPECTER_DIR'] + '/test/data'
        self.psf = load_psf(self.indir + "/psf-spot.fits")

if __name__ == '__main__':
    # unittest.main()           
    s1 = unittest.defaultTestLoader.loadTestsFromTestCase(TestPixPSF)
    s2 = unittest.defaultTestLoader.loadTestsFromTestCase(TestSpotPSF)
    suite = unittest.TestSuite([s1, s2])
    unittest.TextTestRunner(verbosity=2).run(suite)
    # suite = unittest.TestSuite()
    # suite.addTest(TestPixPSF())
    # suite.addTest(TestSpotPSF())
    # suite.run(unittest.TestResult())

