#!/usr/bin/env python

"""
Unit tests for extraction.
"""

import sys
import os
import numpy as np
import scipy.linalg
import unittest
from specter.test import test_data_dir
from specter.psf import load_psf
from specter.extract.ex2d import ex2d, eigen_compose


class TestExtract(unittest.TestCase):
    """
    Test functions within specter.extract
    """
    def setUp(self):
        np.random.seed(0)
        psf = load_psf(test_data_dir() + "/psf-spot.fits")

        nspec = 10
        wmin = min(psf.wavelength(0, y=0), psf.wavelength(nspec-1, y=0))
        ww = psf.wavelength(0, y=np.arange(10,60))
        nwave = len(ww)
        
        phot_shape = (nspec, nwave)
        phot = np.random.uniform(1, 1000, size=phot_shape)
        image_orig = psf.project(ww, phot, verbose=False)
        var = 1.0 + image_orig
        image = image_orig + np.random.normal(scale=np.sqrt(var))
                
        self.phot = phot
        self.image_orig = image_orig
        self.image = image
        self.ivar = 1.0 / var
        self.psf = psf        
        self.ww = ww
        self.nspec = nspec

        # construct a symmetric test matrix
        self.dim = 100
        self.a1 = np.random.uniform(low=0.01, high=100.0, size=(self.dim, self.dim))
        self.a2 = np.random.uniform(low=0.01, high=100.0, size=(self.dim, self.dim))
        self.sym = np.dot(np.transpose(self.a1), self.a1)
        self.sym += np.dot(np.transpose(self.a2), self.a2)

                
    def _test_blat(self):
        from time import time
        specrange = (0, self.nspec)
        waverange = (self.ww[0], self.ww[-1])
        imgvar = 1/self.ivar
        xmin, xmax, ymin, ymax = xyrange = self.psf.xyrange(specrange, waverange)
        
        for i in range(3):
            pix = self.image_orig + np.random.normal(scale=np.sqrt(imgvar))
            d = ex2d(pix, self.ivar, self.psf, specrange, self.ww, full_output=True)
            flux, ivar, R = d['flux'], d['ivar'], d['R']
            rflux = R.dot(self.phot.ravel()).reshape(flux.shape)
            chi = (flux - rflux) * np.sqrt(ivar)
            
            xpix = d['A'].dot(d['xflux'].ravel())
            subpix = pix[ymin:ymax, xmin:xmax].ravel()
            subivar = self.ivar[ymin:ymax, xmin:xmax].ravel()
            
            pixchi = (xpix - subpix) * np.sqrt(subivar)
        
            print i, np.std(chi), np.std(pixchi)
    

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


    def test_noiseless_ex2d(self):
        specrange = (0, self.nspec)
        ivar = np.ones(self.ivar.shape)
        d = ex2d(self.image_orig, ivar, self.psf, specrange, self.ww, full_output=True)

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


    def test_noiseless_ex2d_decorr(self):
        specrange = (0, self.nspec)
        ivar = np.ones(self.ivar.shape)
        d = ex2d(self.image_orig, ivar, self.psf, specrange, self.ww, full_output=True, decorrelate=True)

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


    #- Pull values are wrong.  Why?  Overfitting?
    @unittest.expectedFailure
    def test_ex2d(self):
        specrange = (0, self.nspec)
        d = ex2d(self.image, self.ivar, self.psf, specrange, self.ww, full_output=True)

        #- Pull flux
        R = d['R']
        flux = d['flux']     #- resolution convolved extracted flux
        rphot = R.dot(self.phot.ravel()).reshape(flux.shape)
        pull_flux = (flux - rphot) * np.sqrt(d['ivar'])
        
        #- Pull image
        specrange = (0, self.nspec)
        waverange = (self.ww[0], self.ww[-1])
        xmin, xmax, ymin, ymax = xyrange = self.psf.xyrange(specrange, waverange)
        nx, ny = xmax-xmin, ymax-ymin
        xflux = d['xflux']   #- original extracted flux
        ### ximage = self.psf.project(self.ww, xflux, verbose=False)
        ximage = d['A'].dot(xflux.ravel()).reshape((ny,nx))
        subimg = self.image[ymin:ymax, xmin:xmax]
        subivar = self.ivar[ymin:ymax, xmin:xmax]
        pull_image = ((ximage - subimg) * np.sqrt(subivar))

        print "Known problem: Overfitting may result in small pull value"
        ### print np.std(pull_flux), np.std(pull_image)
        self.assertTrue(np.abs(1-np.std(pull_flux)) < 0.05,
                        msg="pull_flux sigma is %f" % np.std(pull_flux))
        self.assertTrue(np.abs(1-np.std(pull_image)) < 0.05,
                        msg="pull_image sigma is %f" % np.std(pull_image))


    def test_ex2d_subimage(self):
        specrange = (0, self.nspec)
        waverange = self.ww[0], self.ww[-1]
        flux, fluxivar, R = ex2d(self.image, self.ivar, self.psf, specrange, self.ww)

        border = 0
        xmin, xmax, ymin, ymax = self.psf.xyrange(specrange, waverange)
        xmin = max(0, xmin-border)
        xmax = min(self.psf.npix_x, xmax+border)
        ymin = max(0, ymin-border)
        ymax = min(self.psf.npix_y, ymax+border)
        xyrange = (xmin, xmax, ymin, ymax)
        
        subimg = self.image[ymin:ymax, xmin:xmax]
        subivar = self.ivar[ymin:ymax, xmin:xmax]
        subflux, subfluxivar, subR = ex2d(subimg, subivar, self.psf, \
            specrange, self.ww, xyrange=xyrange)

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

        flux, fluxivar, R = ex2d(img, ivar, self.psf, specrange, ww, xyrange=xyrange)
        
        self.assertTrue( np.all(flux == flux) )
        
        
if __name__ == '__main__':
    unittest.main()           
