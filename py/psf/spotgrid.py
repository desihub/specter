#!/usr/bin/env python

"""
SpotGridPSF
"""

import sys
import os
import numpy as N
import fitsio
from specter.psf import PSF
from specter.util import LinearInterp2D, rebin

class SpotGridPSF(PSF):
    def __init__(self, filename):
        """Initialize SpotGridPSF from input file"""
        #- Use PSF class to Load Generic PSF values
        PSF.__init__(self, filename)
        
        #- Load extensions specific to this PSF type
        fx = fitsio.FITS(filename)
        self._spots = fx['SPOTS'].read()
        self._xccd = fx['XCCD'].read()
        self._yccd = fx['YCCD'].read()
        self._fiberpos = fx['FIBERPOS'].read()
        self._spotpos = fx['SPOTPOS'].read()
        self._spotwave = fx['SPOTWAVE'].read()
        
        #- Interpolators
        pp = self._spotpos
        ww = self._spotwave
        self._fspot = LinearInterp2D(pp, ww, self._spots)
        self._fx    = LinearInterp2D(pp, ww, self._xccd)
        self._fy    = LinearInterp2D(pp, ww, self._yccd)
        
        #- Read dimensions from header
        hdr = fx[0].read_header()
        self.CcdPixelSize = hdr['CCDPIXSZ']  #- CCD pixel size in mm
        self.SpotPixelSize = hdr['PIXSIZE']  #- Spot pixel size in mm
        self.FiberSpacing = hdr['DFIBER']   #- center-to-center spacing in mm
        self.GroupSpacing = hdr['DGROUP']   #- center-to-center group gap in mm
        self.FibersPerGroup = hdr['NFIBGRP']
        self.GroupsPerSlit = hdr['NGROUPS']
        self.NumPixX = hdr['NPIX_X']
        self.NumPixY = hdr['NPIX_Y']
        self.nspec = self.FibersPerGroup * self.GroupsPerSlit
        assert(self.nspec == hdr['NSPEC'])
        assert(self.nspec == self._x.shape[0])  #- loaded in PSF.__init__        
        
        fx.close()
        
    def xypix(self, ispec, wavelength):
        """
        Return xslice, yslice, pix for PSF at spectrum ispec, wavelength
        """
        
        #- x,y of spot[0,0], referenced to center of CCD
        p, w = self._fiberpos[ispec], wavelength
        xc = self._fx(p, w)
        yc = self._fy(p, w)
        
        #- Ratio of CCD to Spot pixel sizes
        rpix = int(round(self.CcdPixelSize / self.SpotPixelSize))
        
        #- Calculate offset into CCD pixel
        xoffset = int(xc * rpix) % rpix
        yoffset = int(yc * rpix) % rpix

        #- Place high res spot into grid aligned with CCD pixels
        pix = self._fspot(p, w)
        ny, nx = pix.shape
        A = N.zeros(shape=(pix.shape[0]+rpix, pix.shape[1]+rpix))
        A[yoffset:yoffset+ny, xoffset:xoffset+nx] = pix
        ccdpix = rebin(A, rpix)

        #- TODO: check rounding
        xccd = int(xc - ccdpix.shape[1]/2 + 1)
        yccd = int(yc - ccdpix.shape[0]/2 + 1)
                    
        xx = slice(xccd, xccd+ccdpix.shape[1])
        yy = slice(yccd, yccd+ccdpix.shape[0])
        return xx, yy, ccdpix

        
        
        