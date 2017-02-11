#!/usr/bin/env python

"""
Unit tests for PSF classes.
"""
from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import os
import numpy as np
import unittest
from pkg_resources import resource_filename
from specter.psf import load_psf

class GenericPSFTests(object):
    """
    Wrapper for testing any PSF class
    """

    def wrap_wave_test(self, fn):
        """Test wavelength or loglam"""
        #- Scalar ispec + Unspecified y -> array with npix_y elements
        tmp = fn(0, None)
        self.assertTrue(len(tmp) == self.psf.npix_y)

        #- Scalar ispec and scalar y -> float
        tmp = fn(0, y=0)
        self.assertTrue(isinstance(tmp, float))

        #- Scalar ispec + y array -> wave array
        yy = np.linspace(0, self.psf.npix_y)
        tmp = fn(0, y=yy)
        self.assertTrue(len(tmp) == len(yy))

        #- Unspecified ispec and unspecified y : full wave/loglam grid
        self.assertTrue(fn().shape == (self.psf.nspec, self.psf.npix_y))

        #- ispec >= nspec should raise an error
        with self.assertRaises(IndexError):
            fn(self.psf.nspec)

    #- Test psf.loglam() options
    ### DEPRECATED
    # def test_loglam(self):
    #     self.wrap_wave_test(self.psf.loglam)

    #- Test psf.wavelength() options
    def test_wavelength(self):
        self.wrap_wave_test(self.psf.wavelength)

    #- Make sure that dimensions are set
    def test_dimensions(self):
        self.assertTrue(self.psf.npix_x > 0)
        self.assertTrue(self.psf.npix_y > 0)
        self.assertTrue(self.psf.nspec > 0)

    #- Test xsigma
    def test_xsigma(self):
        yy = (20, self.psf.npix_y//2, self.psf.npix_y-20)
        for ispec in (0, self.psf.nspec//2, self.psf.nspec-1):
            ww = self.psf.wavelength(ispec, y=yy)
            #- Get xsigma for several wavelengths at once
            xsig1 = self.psf.xsigma(ispec, ww)
            self.assertTrue(len(xsig1) == len(ww))
            self.assertTrue(np.min(xsig1) > 0.0)

            #- Call it again to make sure cached results agree
            xsig2 = self.psf.xsigma(ispec, ww)
            self.assertTrue(np.all(xsig1 == xsig2))

        #- Make sure it works for single wavelengths too
        ispec = 0
        ww = self.psf.wavelength(ispec, y=yy)
        xsig1 = self.psf.xsigma(ispec, ww)
        for i in range(len(ww)):
            xsig = self.psf.xsigma(ispec, ww[i])
            self.assertTrue(xsig == xsig1[i])

    #- Test wdisp
    def test_wdisp(self):
        yy = (20, self.psf.npix_y//2, self.psf.npix_y-20)
        for ispec in (0, self.psf.nspec//2, self.psf.nspec-1):
            ww = self.psf.wavelength(ispec, y=yy)
            #- Get wdisp for several wavelengths at once
            xsig1 = self.psf.wdisp(ispec, ww)
            self.assertTrue(len(xsig1) == len(ww))
            self.assertTrue(np.min(xsig1) > 0.0)

            #- Call it again to make sure cached results agree
            xsig2 = self.psf.wdisp(ispec, ww)
            self.assertTrue(np.all(xsig1 == xsig2))

        #- Make sure it works for single wavelengths too
        ispec = 0
        ww = self.psf.wavelength(ispec, y=yy)
        xsig1 = self.psf.wdisp(ispec, ww)
        for i in range(len(ww)):
            xsig = self.psf.wdisp(ispec, ww[i])
            self.assertTrue(xsig == xsig1[i])

    #- Get PSF pixel image at several locations;
    #- Just test that we get a 2D array back
    def test_pix(self):
        ww = self.psf.wavelength()
        wtest = list()
        wtest.append(np.min(ww[:, 0]))
        wtest.append(np.max(ww[:, 0]))
        wtest.append(np.min(ww[:, -1]))
        wtest.append(np.max(ww[:, -1]))
        wtest.append(np.mean(wtest))

        for i in (0, self.psf.nspec//2, self.psf.nspec-1):
            for w in wtest:
                pix = self.psf.pix(i, w)
                self.assertEqual(pix.ndim, 2)

    #- Test that PSF spots are normalized to 1.0
    def test_pix_norm(self):
        psf = self.psf
        yy = np.arange(50, psf.npix_y-50, 50)  #- keep away from the edges
        for i in range(0, psf.nspec, 1+psf.nspec//10):
            ww = psf.wavelength(i, yy)
            for w in ww:
                self.assertAlmostEqual(np.sum(psf.pix(i, w)), 1.0)

    #- Test that PSF spots are positive
    def test_pix_norm(self):
        psf = self.psf
        yy = np.arange(50, psf.npix_y-50, 50)  #- keep away from the edges
        for i in range(0, psf.nspec, 1+psf.nspec//10):
            ww = psf.wavelength(i, yy)
            for w in ww:
                self.assertTrue(np.all(psf.pix(i, w) >= 0), \
                    'PSF goes negative for fiber {} at wavelength {}'.format(i, w))

    #- Get PSF pixel image and where to put it on the CCD
    #- Confirm that image size matches ranges
    def test_xypix(self):
        ww = self.psf.wavelength()
        wtest = list()
        wtest.append(np.min(ww[:, 0]))
        wtest.append(np.max(ww[:, 0]))
        wtest.append(np.min(ww[:, -1]))
        wtest.append(np.max(ww[:, -1]))
        wtest.append(np.mean(wtest))
        wtest.append(np.min(ww)-100)

        for i in (0, self.psf.nspec//2, self.psf.nspec-1):
            for w in wtest:
                xx, yy, pix = self.psf.xypix(self.psf.nspec//2, w)
                shape = (yy.stop-yy.start, xx.stop-xx.start)
                msg = "%s != %s at (i=%d, w=%.1f)" % (str(pix.shape), str(shape), i, w)
                self.assertEqual(pix.shape, shape, msg)

    #- Test psf.xypix() using CCD pixel xmin/xmax, ymin/ymax options
    #- Doesn't test every possibility
    #- TODO: Better tests when walking off edge
    def test_xypix_range(self):
        w = np.mean(self.psf.wavelength())
        i = self.psf.nspec//2
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

    #- Test requests for PSFs off the edge of the requested xyrange
    #- Helper function
    def check_xypix_edges(self, ispec, wavelength, xyrange):
        xx, yy, pix = self.psf.xypix(ispec, wavelength=wavelength, **xyrange)
        nx = xx.stop - xx.start
        ny = yy.stop - yy.start
        msg = 'xx=' + str(xx) + ' yy=' + str(yy) + ' ' + str(pix.shape)
        self.assertEqual(ny, pix.shape[0], msg)
        self.assertEqual(nx, pix.shape[1], msg)

    #- The actual test
    def test_xypix_edges(self):
        psf = self.psf
        #- Pick a range within the CCD
        imin = psf.nspec//2 - 5
        imax = psf.nspec//2 + 5
        wmin = psf.wavelength(imin, y=int(psf.npix_y*0.4))
        wmax = psf.wavelength(imax, y=int(psf.npix_y*0.6))
        wmid = 0.5*(wmin + wmax)

        xmin, xmax, ymin, ymax = psf.xyrange((imin, imax), (wmin, wmax))
        xyrange = dict(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)

        for i in (imin-1, imin, imax, imax+1):
            for w in (wmin-3, wmin, wmid, wmax, wmax+3):
                self.check_xypix_edges(i, w, xyrange)

    #- Test projection of 1D spectrum with 1D wavelength vector
    def test_project11(self):
        ww = self.psf.wavelength(0)[0:10]
        phot = np.random.uniform(0,100,len(ww))
        img = self.psf.project(ww, phot, verbose=False)
        self.assertEqual(img.shape, (self.psf.npix_y, self.psf.npix_x))

    #- Test projection of 2D spectrum with shared 1D wavelength vector
    def test_project12(self):
        ww = self.psf.wavelength(0)[0:10]
        phot = np.random.uniform(0,100,len(ww))
        phot = np.tile(phot, 5).reshape(5, len(ww))
        img = self.psf.project(ww, phot, verbose=False)
        self.assertEqual(img.shape, (self.psf.npix_y, self.psf.npix_x))

    #- Test projection of 2D spectrum with 2D wavelength vector
    def test_project22(self):
        nw = 10
        ww = self.psf.wavelength(0)[0:nw]
        ww = np.tile(ww, 5).reshape(5, nw)
        phot = np.random.uniform(0,100,nw)
        phot = np.tile(phot, 5).reshape(5, nw)
        img = self.psf.project(ww, phot, verbose=False)
        self.assertEqual(img.shape, (self.psf.npix_y, self.psf.npix_x))

    #- Test projection to multiple images at once
    def test_project3(self):
        nwave = 10
        nspec = 5
        nimg = 3
        ww = self.psf.wavelength(0)[0:nwave]
        ww = np.tile(ww, nspec).reshape(nspec, nwave)
        phot = np.random.uniform(0,100,size=(nimg, nspec, nwave))
        img = self.psf.project(ww, phot, verbose=False)
        self.assertEqual(img.shape, (nimg, self.psf.npix_y, self.psf.npix_x))
        img0 = self.psf.project(ww, phot[0], verbose=False)
        img1 = self.psf.project(ww, phot[1], verbose=False)
        self.assertTrue(np.all(img0==img[0]))
        self.assertTrue(np.all(img1==img[1]))

        #- Confirm that it also works with 1D wave and 3D phot
        imgx = self.psf.project(ww[0], phot, verbose=False)
        self.assertTrue(np.all(img==imgx))

    #- Test projection starting at specmin != 0
    def test_project_specmin(self):
        ww = self.psf.wavelength(0)[0:10]
        phot = np.random.uniform(0,100,len(ww))
        img = self.psf.project(ww, phot, specmin=1, verbose=False)
        self.assertEqual(img.shape, (self.psf.npix_y, self.psf.npix_x))

        #- specmin >= nspec should raise an error
        with self.assertRaises(ValueError):
            i = self.psf.nspec
            img = self.psf.project(ww, phot, specmin=i, verbose=False)

    #- Test projecting to a subgrid of CCD pixels
    def test_project_xyrange(self):
        nspec = 5
        nw = 10
        ww = self.psf.wavelength(0)[0:nw]
        phot = np.random.uniform(0,100,nw)               #- 1D
        phot = np.tile(phot, nspec).reshape(nspec, nw)   #- 2D

        spec_range = (0, nspec)

        xyrange = xmin,xmax,ymin,ymax = self.psf.xyrange(spec_range, ww)
        img = self.psf.project(ww, phot, verbose=False)
        subimg = self.psf.project(ww, phot, xyrange=xyrange, verbose=False)

        #- Does subimage match corresponding range for full image?
        self.assertTrue(np.all(subimg == img[ymin:ymax, xmin:xmax]))

        #- Clear subimage region and test that everything is 0
        img[ymin:ymax, xmin:xmax] = 0.0
        self.assertTrue(np.all(img == 0.0))

    #- Test projection with an xyrange that is smaller than wavelength range
    def test_project_small_xyrange(self):
        #- Find the xyrange for a small range of wavelengths
        nspec = 5
        nw = 5
        ww = self.psf.wavelength(0)[1000:1000+nw]
        spec_range = (0, nspec)
        xyrange = xmin,xmax,ymin,ymax = self.psf.xyrange(spec_range, ww)

        #- Now create a spectrum with larger wavelength coverage
        ww = self.psf.wavelength(0)[1000-5*nw:1000+5*nw]
        nw = len(ww)
        phot = np.random.uniform(0,100,nw)               #- 1D
        phot = np.tile(phot, nspec).reshape(nspec, nw)   #- 2D

        subimg = self.psf.project(ww, phot, xyrange=xyrange, verbose=False)

    #- Test the projection matrix gives same answer as psf.project()
    def test_projection_matrix(self):
        nspec = 5
        nw = 20
        w_edge = 10  #- avoid edge effects; test that separately
        phot = np.random.uniform(100,1000, size=(nspec, nw))
        for specmin in (0, self.psf.nspec//2, self.psf.nspec-nspec-1):
            specrange = (specmin, specmin+nspec)
            wspec = self.psf.wavelength(specmin, [0, self.psf.npix_y])
            for wmin in (wspec[0]+w_edge, 0.5*(wspec[0]+wspec[1]), wspec[1]-nw-w_edge):
                ww = np.arange(wmin, wmin+nw)
                waverange = (ww[0], ww[-1])
                xmin, xmax, ymin, ymax = xyrange = self.psf.xyrange(specrange, waverange)
                nx = xmax-xmin
                ny = ymax-ymin

                img1 = self.psf.project(ww, phot, xyrange=xyrange, \
                                        specmin=specmin, verbose=False)

                A = self.psf.projection_matrix(specrange, ww, xyrange)
                img2 = A.dot(phot.ravel()).reshape((ny, nx))

                self.assertTrue(np.all(img1==img2))

    #- Test xyrange with scalar vs. tuple spec_range
    def test_xyrange_ispec(self):
        ispec = 0
        ww = self.psf.wavelength(ispec, y=np.arange(0, 10))
        xyr1 = self.psf.xyrange(ispec, ww)
        xyr2 = self.psf.xyrange((ispec, ispec+1), ww)
        self.assertTrue(xyr1 == xyr2)

    #- Test xyrange on edges of CCD
    def test_xyrange_edges(self):
        psf = self.psf
        y = np.arange(-1, psf.npix_y, 0.2)
        for i in range(0, psf.nspec, psf.nspec//10):
            w = psf.wavelength(i, y)
            xmin, xmax, ymin, ymax = psf.xyrange(i, w)
            self.assertTrue(xmin <= xmax)
            self.assertTrue(ymin <= ymax)

    #- Test projection matrix with scalar vs. tuple spec_range
    def test_projmat_ispec(self):
        ispec = 0
        ww = self.psf.wavelength(ispec, y=np.arange(0, 10))
        xyrange = self.psf.xyrange(ispec, ww)
        A = self.psf.projection_matrix(ispec, ww, xyrange)
        B = self.psf.projection_matrix((ispec, ispec+1), ww, xyrange)

        self.assertTrue(A.shape == B.shape)
        self.assertTrue(np.all(A.data == B.data))

    #- Test shift of PSF xy solution
    @unittest.expectedFailure
    def test_shift_xy(self):
        dx, dy = 0.1, 0.2
        x0 = self.psf.x(0).copy()
        y0 = self.psf.y(0).copy()
        self.psf.shift_xy(dx, dy)
        self.assertTrue(np.all(self.psf.x(0) == x0+dx))
        self.assertTrue(np.all(self.psf.y(0) == y0+dy))

    #- Test multiple options for getting x centroid
    def test_x(self):
        #- Grid of x positions
        x = self.psf.x()
        self.assertEqual(x.ndim, 2)
        self.assertEqual(x.shape[0], self.psf.nspec)
        self.assertEqual(x.shape[1], self.psf.npix_y)

        #- x for fiber 0
        x = self.psf.x(0)
        self.assertEqual(x.ndim, 1)
        self.assertEqual(len(x), self.psf.npix_y)

        #- x for fiber 0 at a specific wavelength
        w = np.mean(self.psf.wavelength(0))
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

    #- Test multiple options for getting y centroid
    def test_y(self):
        #- Grid of y positions
        # y = self.psf.y()
        # self.assertEqual(y.ndim, 2)
        # self.assertEqual(y.shape[0], self.psf.nspec)
        # self.assertEqual(y.shape[1], self.psf.npix_y)

        #- y for fiber 0
        # y = self.psf.y(0)
        # self.assertEqual(y.ndim, 1)
        # self.assertEqual(len(y), self.psf.npix_y)

        #- y for fiber 0 at a specific wavelength
        w = np.mean(self.psf.wavelength(0))
        y = self.psf.y(0, w)
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

    #- Ensure that pix requests outside of wavelength range are blank
    def test_waverange_pix(self):
        for ispec in (0, 1, self.psf.nspec//2, self.psf.nspec-1):
            xx, yy, pix = self.psf.xypix(ispec, self.psf.wmin-1)
            self.assertTrue(xx.start == xx.stop == 0)
            self.assertTrue(yy.start == yy.stop == 0)
            self.assertTrue(pix.shape == (0,0))

            xx, yy, pix = self.psf.xypix(ispec, self.psf.wmax+1)
            self.assertTrue(xx.start == xx.stop == 0)
            self.assertTrue(yy.start == yy.stop == self.psf.npix_y)
            self.assertTrue(pix.shape == (0,0))


    #- Test getting x and y at the same time
    def test_xy(self):
        w = self.psf.wavelength(0)
        x = self.psf.x(0, w)
        y = self.psf.y(0, w)
        xy = self.psf.xy(0, w)
        self.assertTrue(np.all(xy[0] == x))
        self.assertTrue(np.all(xy[1] == y))

    #- sanity check on wmin, wmax, wmin_all, wmax_all
    def test_waveminmax(self):
        psf = self.psf
        self.assertLessEqual(psf.wmin, psf.wmin_all)
        self.assertGreaterEqual(psf.wmax, psf.wmax_all)
        self.assertLess(psf.wmin, psf.wmax)

    #- Test getting x and y and wavelength at the same time
    ### REMOVED
    # def test_xyw(self):
    #     x = self.psf.x(0)
    #     y = self.psf.y(0)
    #     w = self.psf.wavelength(0)
    #     xyw = self.psf.xyw(0)
    #     self.assertTrue(np.all(xyw[0] == x))
    #     self.assertTrue(np.all(xyw[1] == y))
    #     self.assertTrue(np.all(xyw[2] == w))

    #- Test getting xy pixel range where spectra would project
    #+ TODO: Test case where input wavelength grid is sub-pixel
    def test_xyrange(self):
        ww = self.psf.wavelength(0)[0:10]
        spec_range = [0,10]
        xmin,xmax,ymin,ymax = self.psf.xyrange(spec_range, ww)

        #- Test all wavelengths for first and last spectrum
        for i in range(spec_range[0], spec_range[-1]):
            for w in ww:
                xx, yy, pix = self.psf.xypix(i, w)
                if pix.shape == (0,0):  #- off edge of CCD
                    continue

                self.assertGreaterEqual(xx.start, xmin)
                self.assertLessEqual(xx.stop, xmax)
                self.assertGreaterEqual(yy.start, ymin)
                self.assertLessEqual(yy.stop, ymax)

        #- Test all spectra for min and max wavelengths
        for i in range(spec_range[0], spec_range[-1]):
            for w in (ww[0], ww[-1]):
                xx, yy, pix = self.psf.xypix(i, w)
                if pix.shape == (0,0):  #- off edge of CCD
                    continue

                self.assertGreaterEqual(xx.start, xmin)
                self.assertLessEqual(xx.stop, xmax)
                self.assertGreaterEqual(yy.start, ymin)
                self.assertLessEqual(yy.stop, ymax)

        #- Test that extending off the chip is still ok
        wmin = self.psf.wmin - 10
        wmax = self.psf.wmax + 10
        xmin, xmax, ymin, ymax = self.psf.xyrange(spec_range, (wmin, wmax))
        self.assertEqual(ymin, 0)
        self.assertEqual(ymax, self.psf.npix_y)
        self.assertLess(xmin, xmax)
        self.assertGreater(xmax, 0)

#- Test Pixellated PSF format
class TestPixPSF(GenericPSFTests,unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.psf = load_psf(resource_filename("specter.test", "t/psf-pix.fits"))

#- Test SpotGrid PSF format
class TestSpotPSF(GenericPSFTests,unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.psf = load_psf(resource_filename("specter.test", "t/psf-spot.fits"))

#- Test MonoSpot PSF format
class TestMonoSpotPSF(GenericPSFTests,unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.psf = load_psf(resource_filename("specter.test", "t/psf-monospot.fits"))

#- Test GaussHermitePSF format
class TestGaussHermitePSF(GenericPSFTests,unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.psf = load_psf(resource_filename("specter.test", "t/psf-gausshermite.fits"))

    def test_spot_center(self):
        import scipy.ndimage
        for y in np.linspace(14,16,9):
            w = self.psf.wavelength(0, y)
            x, y = self.psf.xy(0, w)
            xx, yy, pix = self.psf.xypix(0, w)

            y0, x0 = scipy.ndimage.center_of_mass(pix)
            dx = x0 - pix.shape[1]//2
            dy = y0 - pix.shape[0]//2

            self.assertLessEqual(np.abs(dx), 0.55)
            self.assertLessEqual(np.abs(dy), 0.55)
            self.assertLessEqual(np.abs(xx.start+x0-x), 0.1)
            self.assertLessEqual(np.abs(yy.start+y0-y), 0.1)

#- Test GaussHermite2PSF format
class TestGaussHermite2PSF(GenericPSFTests,unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.psf = load_psf(resource_filename("specter.test", "t/psf-gausshermite2.fits"))

if __name__ == '__main__':

    import warnings
    warnings.simplefilter('error')

    # unittest.main()
    s1 = unittest.defaultTestLoader.loadTestsFromTestCase(TestPixPSF)
    s2 = unittest.defaultTestLoader.loadTestsFromTestCase(TestSpotPSF)
    s3 = unittest.defaultTestLoader.loadTestsFromTestCase(TestMonoSpotPSF)
    s4 = unittest.defaultTestLoader.loadTestsFromTestCase(TestGaussHermitePSF)
    s5 = unittest.defaultTestLoader.loadTestsFromTestCase(TestGaussHermite2PSF)

    unittest.TextTestRunner(verbosity=2).run(unittest.TestSuite([s1, s2, s3, s4, s5]))
    # unittest.TextTestRunner(verbosity=2).run(unittest.TestSuite([s2,]))
