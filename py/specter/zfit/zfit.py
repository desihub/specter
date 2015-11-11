#!/usr/bin/env python

"""
IN PROGRESS; Fit for redshift given a template
"""

import sys
import os
import numpy as np
from scipy.ndimage.filters import gaussian_filter
from specter.util import resample

def zchi2(spectra, wt, template, z):
    wz = wt*(1+z)
    chi2 = 0.0
    for (wave, flux, ivar, psf, ispec) in spectra:
        #- Wavelength at each CCD row
        wpsf = psf.wavelength(ispec, y=np.arange(psf.npix_y))

        #- Template pixel size in Angstroms at those wavelengths
        dw = np.interp(wpsf, wz, np.gradient(wz))        
        
        #- wavelength dispersion in units of template pixels
        wdisp = np.median( psf.wdisp(wpsf) / dw )
        
        #- Convolve with resolution and then sample with input spectrum binning
        model = resample(wz, gaussian_filter(template, wdisp_pix), xnew=wave)

        chi2 += np.sum( (flux-model)**2 * ivar )
        
    return chi2

def zfit(spectra, wt, template, zmin=0.0, zmax=2.0, dlogz=0.001):
    """
    spectra: array of spectra, each of which is a (wave, flux, ivar, psf, ispec) tuple
    wt : wavelength array of template spectrum
    template : template flux
    """
    
    zz = 10**np.arange(np.log10(1.0+zmin), np.log10(1.0+zmax), dlogz)
    chi2 = np.zeros(len(zz))
    for i in range(zz):
        print z[i], i, len(zz)
        chi2[i] = zchi2(spectra, wt, template, zz[i])
        
    return zz, chi2
    
    
    
    
    