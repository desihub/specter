#!/usr/bin/env python

"""
Unit tests for PSF classes.
"""

import sys
import os
import numpy as N
import unittest
from specter.test import test_data_dir
from specter.psf import load_psf
from specter.extract.ex2d import ex2d


class TestExtract(unittest.TestCase):
    """
    Test functions within specter.util
    """
    def setUp(self):
        N.random.seed(0)
        psf = load_psf(test_data_dir() + "/psf-spot.fits")

        nspec = 10
        wmin = min(psf.wavelength(0, y=0), psf.wavelength(nspec-1, y=0))
        ww = N.arange(wmin+10, wmin+60)
        nwave = len(ww)
        
        phot_shape = (nspec, nwave)
        phot = N.random.uniform(100,1000, size=phot_shape)
        image_orig = psf.project(phot, ww, verbose=False)
        var = 1.0 + N.abs(image_orig)
        image = image_orig + N.random.normal(scale=N.sqrt(var))
                
        self.phot = phot
        self.image_orig = image_orig
        self.image = image
        self.ivar = 1.0 / var
        self.psf = psf        
        self.ww = ww
        self.nspec = nspec
                
    
    def test_noiseless_ex2d(self):
        specrange = (0, self.nspec)
        ivar = N.ones(self.ivar.shape)
        d = ex2d(self.image_orig, ivar, self.psf, specrange, self.ww, full_output=True)

        R = d['R']
        flux = d['flux']     #- resolution convolved extracted flux
        xflux = d['xflux']   #- original extracted flux
        
        #- Resolution convolved input photons (flux)
        rphot = R.dot(self.phot.ravel()).reshape(flux.shape)
        
        #- extracted flux projected back to image
        ximg = self.psf.project(xflux, self.ww, verbose=False)
        
        #- Compare inputs to outputs
        dflux = flux - rphot
        dxflux = xflux - self.phot
        dximg = ximg - self.image_orig
        self.assertTrue( N.max(N.abs(dflux)) < 1e-9 )
        self.assertTrue( N.max(N.abs(dxflux)) < 1e-9 )
        self.assertTrue( N.max(N.abs(dximg)) < 1e-9 )                

    def test_ex2d(self):
        specrange = (0, self.nspec)
        d = ex2d(self.image, self.ivar, self.psf, specrange, self.ww, full_output=True)

        #- Pull flux
        R = d['R']
        flux = d['flux']     #- resolution convolved extracted flux
        rphot = R.dot(self.phot.ravel()).reshape(flux.shape)
        pull_flux = (flux - rphot) * N.sqrt(d['ivar'])
        
        #- Pull image
        specrange = (0, self.nspec)
        waverange = (self.ww[0], self.ww[-1])
        xmin, xmax, ymin, ymax = xyrange = self.psf.xyrange(specrange, waverange)
        xflux = d['xflux']   #- original extracted flux
        ximage = self.psf.project(xflux, self.ww, verbose=False)
        pull_image = ((ximage - self.image) * N.sqrt(self.ivar))[ymin:ymax, xmin:xmax]
        
        
        self.assertTrue(N.abs(1-N.std(pull_flux)) < 0.03,
                        msg="pull_flux sigma is %f" % N.std(pull_flux))
        self.assertTrue(N.abs(1-N.std(pull_image)) < 0.03,
                        msg="pull_image sigma is %f" % N.std(pull_image))
        
    def test_ex2d_subimage(self):
        specrange = (0, self.nspec)
        waverange = self.ww[0], self.ww[-1]
        flux, ivar, R = ex2d(self.image, self.ivar, self.psf, specrange, self.ww)

        xmin, xmax, ymin, ymax = self.psf.xyrange(specrange, waverange)
        xmin = max(0, xmin-10)
        xmax = min(self.psf.npix_x, xmax+10)
        ymin = max(0, ymin-10)
        ymax = min(self.psf.npix_y, ymax+10)
        xyrange = (xmin, xmax, ymin, ymax)
        
        subimg = self.image[ymin:ymax, xmin:xmax]
        subivar = self.ivar[ymin:ymax, xmin:xmax]
        subflux, subivar, subR = ex2d(subimg, subivar, self.psf, \
            specrange, self.ww, xyrange=xyrange)

        self.assertTrue( N.all(subflux == flux) )
        self.assertTrue( N.all(subivar == ivar) )
        self.assertTrue( N.all(subR == R) )
        
if __name__ == '__main__':
    unittest.main()           
