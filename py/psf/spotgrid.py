#!/usr/bin/env python
"""
SpotGridPSF - Linear interpolate hi-res sampled spots to model PSF

Stephen Bailey
Fall 2012
"""

import sys
import os
import numpy as N
import fitsio
from specter.psf import PSF
from specter.util import LinearInterp2D, rebin

class SpotGridPSF(PSF):
    """
    Model PSF with a linear interpolation of high resolution sampled spots
    """
    def __init__(self, filename):
        """
        Initialize SpotGridPSF from input file
        
        See specter.psf.PSF for futher details
        """
        #- Use PSF class to Load Generic PSF values (x, y, wavelength, ...)
        PSF.__init__(self, filename)
        
        #- Load extensions specific to this PSF type
        fx = fitsio.FITS(filename)
        self._spots = fx['SPOTS'].read()  #- PSF spots
        self._spotx  = fx['SPOTX'].read()   #- X location of spots
        self._spoty  = fx['SPOTY'].read()   #- Y location of spots
        self._fiberpos = fx['FIBERPOS'].read()  #- Location of fibers on slit
        self._spotpos = fx['SPOTPOS'].read()    #- Slit loc of sampled spots
        self._spotwave = fx['SPOTWAVE'].read()  #- Wavelengths of spots
        
        #- 2D linerar interpolators
        pp = self._spotpos
        ww = self._spotwave
        self._fspot = LinearInterp2D(pp, ww, self._spots)
        self._fx    = LinearInterp2D(pp, ww, self._spotx)
        self._fy    = LinearInterp2D(pp, ww, self._spoty)
        
        #- Read spot vs. CCD pixel scales from header
        hdr = fx[0].read_header()
        self.CcdPixelSize = hdr['CCDPIXSZ']  #- CCD pixel size in mm
        self.SpotPixelSize = hdr['PIXSIZE']  #- Spot pixel size in mm
        
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

        #- TODO: Could do sub-pixel sinc shifting, but that is slow

        #- Find where the [0,0] pixel goes on the CCD 
        xccd = int(xc - ccdpix.shape[1]/2 + 1)
        yccd = int(yc - ccdpix.shape[0]/2 + 1)
        
        #- Check if completely off the edge in any direction
        if (xccd > self.npix_x) or (xccd+ccdpix.shape[1] < 0) or \
           (yccd > self.npix_y) or (yccd+ccdpix.shape[0] < 0):
            return slice(0,0), slice(0,0), N.zeros(0)
            
        #- Check if partially off edge
        if xccd < 0:
            ccdpix = ccdpix[:, -xccd:]
            xccd = 0
        elif xccd + ccdpix.shape[1] > self.npix_x:
            dx = xccd + ccdpix.shape[1] - self.npix_x
            ccdpix = ccdpix[:, -dx]

        if yccd < 0:
            ccdpix = ccdpix[-yccd:, ]
            yccd = 0
        elif yccd + ccdpix.shape[0] > self.npix_y:
            dy = yccd + ccdpix.shape[0] - self.npix_y
            ccdpix = ccdpix[-dy, :]
        
        xx = slice(xccd, xccd+ccdpix.shape[1])
        yy = slice(yccd, yccd+ccdpix.shape[0])
        
        return xx, yy, ccdpix

        
        
        