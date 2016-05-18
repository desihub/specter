#!/usr/bin/env python

"""
Unit tests for extraction.
"""

import sys
import os
import numpy as np
import scipy.linalg
import unittest
from pkg_resources import resource_filename
from ..psf import load_psf
from ..extract.ex2d import ex2d, ex2d_patch, eigen_compose
from ..extract.ex1d import ex1d

class TestExtract(unittest.TestCase):
    """
    Test functions within specter.extract
    """
    @classmethod
    def setUpClass(cls):
        cls.psf = load_psf(resource_filename("specter.test", "t/psf-spot.fits"))

        np.random.seed(0)
        nspec = 10
        wmin = min(cls.psf.wavelength(0, y=0), cls.psf.wavelength(nspec-1, y=0))
        ### ww = psf.wavelength(0, y=np.arange(10,60))
        wmin, wmax = cls.psf.wavelength(0, y=(10,50))
        ww = np.arange(wmin, wmax, 0.5)
        nwave = len(ww)

        phot_shape = (nspec, nwave)
        phot = np.random.uniform(1, 1000, size=phot_shape)
        image_orig = cls.psf.project(ww, phot, verbose=False)
        var = 1.0 + image_orig
        image = image_orig + np.random.normal(scale=np.sqrt(var))

        cls.phot = phot
        cls.image_orig = image_orig
        cls.image = image
        cls.ivar = 1.0 / var
        cls.ww = ww
        cls.nspec = nspec

        # construct a symmetric test matrix
        cls.dim = 100
        cls.a1 = np.random.uniform(low=0.01, high=100.0, size=(cls.dim, cls.dim))
        cls.a2 = np.random.uniform(low=0.01, high=100.0, size=(cls.dim, cls.dim))
        cls.sym = np.dot(np.transpose(cls.a1), cls.a1)
        cls.sym += np.dot(np.transpose(cls.a2), cls.a2)


    def test_ex2d_chi2(self):
        from time import time
        specrange = (0, self.nspec)
        waverange = (self.ww[0], self.ww[-1])
        imgvar = 1/self.ivar
        xmin, xmax, ymin, ymax = xyrange = self.psf.xyrange(specrange, waverange)

        for i in range(3):
            pix = self.image_orig + np.random.normal(scale=np.sqrt(imgvar))
            d = ex2d_patch(pix, self.ivar, self.psf, 0, self.nspec, self.ww, full_output=True)
            flux, ivar, R = d['flux'], d['ivar'], d['R']
            rflux = R.dot(self.phot.ravel()).reshape(flux.shape)
            chi = (flux - rflux) * np.sqrt(ivar)
            #- loose test, just checking for catastrophic failures
            self.assertLess(abs(1-np.std(chi)), 0.10)

    def test_eigen_compose(self):
        w, v = scipy.linalg.eigh(self.sym)
        check = eigen_compose(w, v)
        np.testing.assert_almost_equal(check, self.sym, decimal=5)
        check = eigen_compose(w, v, sqr=True)
        check = eigen_compose(w, v, invert=True)
        check = eigen_compose(w, v, invert=True, sqr=True)

        # check reconstruction
        w_inv, v_inv = scipy.linalg.eigh(check)
        comp_w = np.multiply(w_inv, w_inv)
        comp = eigen_compose(comp_w, v_inv, invert=True)
        np.testing.assert_almost_equal(comp, self.sym, decimal=3)

    def test_ex1d(self):
        specrange = (0, self.nspec)
        mask = np.zeros(self.image.shape, dtype=int)
        ny = 20
        flux, ivar = ex1d(self.image, mask, self.psf, yrange=[0,ny],
            readnoise=1.0, specrange=specrange)
        self.assertEqual(flux.shape, ivar.shape)
        self.assertEqual(flux.shape[0], self.nspec)
        self.assertEqual(flux.shape[1], ny)

    def test_ex1d_model(self):
        specrange = (2, self.nspec)
        mask = np.zeros(self.image.shape, dtype=int)
        ny = 20
        flux, ivar, model = ex1d(self.image, mask, self.psf, yrange=[ny,2*ny],
            readnoise=1.0, specrange=specrange, model=True)
        self.assertEqual(flux.shape, ivar.shape)
        self.assertEqual(flux.shape[0], self.nspec-2)
        self.assertEqual(flux.shape[1], ny)
        self.assertEqual(self.image.shape, model.shape)

    def test_noiseless_ex2d(self):
        specrange = (0, self.nspec)
        ivar = np.ones(self.ivar.shape)
        d = ex2d_patch(self.image_orig, ivar, self.psf, 0, self.nspec, self.ww, full_output=True, ndecorr=True)

        R = d['R']
        flux = d['flux']     #- resolution convolved extracted flux
        xflux = d['xflux']   #- original extracted flux

        #- Resolution convolved input photons (flux)
        rphot = R.dot(self.phot.ravel()).reshape(flux.shape)

        #- extracted flux projected back to image
        ximg = self.psf.project(self.ww, xflux, verbose=False)

        #- Compare inputs to outputs
        bias = (flux - rphot)/rphot
        dximg = ximg - self.image_orig

        self.assertTrue( np.max(np.abs(bias)) < 1e-9 )
        self.assertTrue( np.max(np.abs(dximg)) < 1e-6 )


    def test_noiseless_ex2d_sigdecorr(self):
        specrange = (0, self.nspec)
        ivar = np.ones(self.ivar.shape)
        d = ex2d_patch(self.image_orig, ivar, self.psf, 0, self.nspec, self.ww, full_output=True)

        R = d['R']
        flux = d['flux']     #- resolution convolved extracted flux
        xflux = d['xflux']   #- original extracted flux

        #- Resolution convolved input photons (flux)
        rphot = R.dot(self.phot.ravel()).reshape(flux.shape)

        #- extracted flux projected back to image
        ximg = self.psf.project(self.ww, xflux, verbose=False)

        #- Compare inputs to outputs
        bias = (flux - rphot)/rphot
        dximg = ximg - self.image_orig

        self.assertTrue( np.max(np.abs(bias)) < 1e-9 )
        self.assertTrue( np.max(np.abs(dximg)) < 1e-6 )

    def test_ex2d(self):
        flux, ivar, Rdata = ex2d(self.image, self.ivar, self.psf, 0, self.nspec,
            self.ww, wavesize=len(self.ww)//5)
        self.assertEqual(flux.shape, (self.nspec, len(self.ww)))
        self.assertEqual(ivar.shape, (self.nspec, len(self.ww)))

        self.assertEqual(Rdata.ndim, 3)
        self.assertEqual(Rdata.shape[0], self.nspec)
        self.assertEqual(Rdata.shape[2], len(self.ww))
        self.assertGreater(Rdata.shape[1], 15)

    #- Pull values are wrong.  Why?  Overfitting?
    def test_ex2d_patch(self):
        d = ex2d_patch(self.image, self.ivar, self.psf, 0, self.nspec, self.ww, full_output=True)

        self.assertEqual(d['flux'].shape, (self.nspec, len(self.ww)))
        self.assertEqual(d['ivar'].shape, (self.nspec, len(self.ww)))

        ntot = len(self.ww) * self.nspec
        self.assertEqual(d['R'].shape, (ntot, ntot))

        #- Pull flux
        # R = d['R']
        # flux = d['flux']     #- resolution convolved extracted flux
        # rphot = R.dot(self.phot.ravel()).reshape(flux.shape)
        # pull_flux = (flux - rphot) * np.sqrt(d['ivar'])
        #
        # #- Pull image
        # specrange = (0, self.nspec)
        # waverange = (self.ww[0], self.ww[-1])
        # xmin, xmax, ymin, ymax = xyrange = self.psf.xyrange(specrange, waverange)
        # nx, ny = xmax-xmin, ymax-ymin
        # xflux = d['xflux']   #- original extracted flux
        # ### ximage = self.psf.project(self.ww, xflux, verbose=False)
        # ximage = d['A'].dot(xflux.ravel()).reshape((ny,nx))
        # subimg = self.image[ymin:ymax, xmin:xmax]
        # subivar = self.ivar[ymin:ymax, xmin:xmax]
        # pull_image = ((ximage - subimg) * np.sqrt(subivar))

        # print "Known problem: Overfitting may result in small pull value"
        # ### print np.std(pull_flux), np.std(pull_image)
        # self.assertTrue(np.abs(1-np.std(pull_flux)) < 0.05,
        #                 msg="pull_flux sigma is %f" % np.std(pull_flux))
        # self.assertTrue(np.abs(1-np.std(pull_image)) < 0.05,
        #                 msg="pull_image sigma is %f" % np.std(pull_image))


    def test_ex2d_subimage(self):
        specrange = (0, self.nspec)
        waverange = self.ww[0], self.ww[-1]
        flux, fluxivar, R = ex2d_patch(self.image, self.ivar, self.psf, 0, self.nspec, self.ww)

        border = 0
        xmin, xmax, ymin, ymax = self.psf.xyrange(specrange, waverange)
        xmin = max(0, xmin-border)
        xmax = min(self.psf.npix_x, xmax+border)
        ymin = max(0, ymin-border)
        ymax = min(self.psf.npix_y, ymax+border)
        xyrange = (xmin, xmax, ymin, ymax)

        subimg = self.image[ymin:ymax, xmin:xmax]
        subivar = self.ivar[ymin:ymax, xmin:xmax]
        subflux, subfluxivar, subR = ex2d_patch(subimg, subivar, self.psf, \
            0, self.nspec, self.ww, xyrange=xyrange)

        #- These arrays pass np.allclose, but sometimes fail np.all on edison.
        #- They should be absolutely identical.  Leaving this as failing.
        self.assertTrue( np.all(subflux == flux) )
        self.assertTrue( np.all(subfluxivar == fluxivar) )
        self.assertTrue( np.all(subR == R) )


    def test_wave_off_image(self):
        ww = self.psf.wmin - 5 + np.arange(10)
        nspec = 2
        specrange = [0,nspec]
        xyrange = self.psf.xyrange(specrange, ww)

        phot = np.ones([nspec,len(ww)])

        img = self.psf.project(ww, phot, xyrange=xyrange)
        ivar = np.ones(img.shape)

        flux, fluxivar, R = ex2d_patch(img, ivar, self.psf, 0, self.nspec, ww, xyrange=xyrange)

        self.assertTrue( np.all(flux == flux) )


if __name__ == '__main__':
    unittest.main()
