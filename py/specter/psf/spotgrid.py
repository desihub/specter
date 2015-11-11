#!/usr/bin/env python
"""
SpotGridPSF - Linear interpolate hi-res sampled spots to model PSF

Stephen Bailey
Fall 2012
"""

import sys
import os
import numpy as np
from astropy.io import fits
from specter.psf import PSF
from specter.util import LinearInterp2D, rebin_image, sincshift
import scipy.interpolate
import sys # for debug
import pylab # for debug

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
        fx = fits.open(filename)
        self._spots = fx['SPOTS'].data        #- PSF spots
        # self._spotx  = fx['SPOTX'].data     #- X location of spots
        # self._spoty  = fx['SPOTY'].data     #- Y location of spots
        self._fiberpos = fx['FIBERPOS'].data  #- Location of fibers on slit
        self._spotpos = fx['SPOTPOS'].data    #- Slit loc of sampled spots
        self._spotwave = fx['SPOTWAVE'].data  #- Wavelengths of spots
        
        #- 2D linerar interpolators
        pp = self._spotpos
        ww = self._spotwave
        self._fspot = LinearInterp2D(pp, ww, self._spots)
        # self._fx    = LinearInterp2D(pp, ww, self._spotx)
        # self._fy    = LinearInterp2D(pp, ww, self._spoty)
        
        #- Read spot vs. CCD pixel scales from header
        hdr = fx[0].header
        self.CcdPixelSize = hdr['CCDPIXSZ']  #- CCD pixel size in mm
        self.SpotPixelSize = hdr['PIXSIZE']  #- Spot pixel size in mm
        
        fx.close()
        
    def _xypix(self, ispec, wavelength):
        """
        Return xslice, yslice, pix for PSF at spectrum ispec, wavelength
        """
        #return self._xypix_sincshift(ispec, wavelength)
        return self._xypix_histo(ispec, wavelength)

    def bilinear_interpolate(self,im, x, y):
        x = N.asarray(x)
        y = N.asarray(y)

        x0 = N.floor(x).astype(int)
        x1 = x0 + 1
        y0 = N.floor(y).astype(int)
        y1 = y0 + 1

        x0 = N.clip(x0, 0, im.shape[1]-1);
        x1 = N.clip(x1, 0, im.shape[1]-1);
        y0 = N.clip(y0, 0, im.shape[0]-1);
        y1 = N.clip(y1, 0, im.shape[0]-1);

        Ia = im[ y0, x0 ]
        Ib = im[ y1, x0 ]
        Ic = im[ y0, x1 ]
        Id = im[ y1, x1 ]
        
        wa = (x1-x) * (y1-y)
        wb = (x1-x) * (y-y0)
        wc = (x-x0) * (y1-y)
        wd = (x-x0) * (y-y0)

        return wa*Ia + wb*Ib + wc*Ic + wd*Id
    
    def _xypix_histo(self, ispec, wavelength):
        """
        Return xslice, yslice, pix for PSF at spectrum ispec, wavelength
        """
        print "in _xypix_histo"
        
        #- x,y of spot on CCD
        p, w = self._fiberpos[ispec], wavelength
        x_ccd_center, y_ccd_center = self.xy(ispec, wavelength)
        print "x_ccd_center, y_ccd_center = ",x_ccd_center, y_ccd_center
        
        #- Ratio of CCD to Spot pixel sizes
        rebin = (self.CcdPixelSize / self.SpotPixelSize)
        if rebin-int(rebin)>1e-6 :
            print "ERROR I CANT DO A SIMPLE REBINNING"
            sys.exit(12)
        rebin=int(rebin)
        pix_spot_values=self._fspot(p, w)
        nx_spot=pix_spot_values.shape[1]
        ny_spot=pix_spot_values.shape[0]
        print nx_spot,ny_spot
        dx=(x_ccd_center*rebin-int(x_ccd_center*rebin))/rebin
        dy=(y_ccd_center*rebin-int(y_ccd_center*rebin))/rebin
        tmp_x=N.tile(N.arange(nx_spot)-dx,(ny_spot,1))
        tmp_y=N.tile(N.arange(ny_spot)-dx,(nx_spot,1)).T # is it ok ?
        
        resampled_pix_spot_values=self.bilinear_interpolate(pix_spot_values,tmp_x,tmp_y)

        nx_ccd=nx_spot/rebin
        ny_ccd=ny_spot/rebin
        
        print pix_spot_values.shape
        print resampled_pix_spot_values.shape
        
        ccd_pix_spot_values=resampled_pix_spot_values.reshape(ny_spot,nx_ccd,rebin).sum(2)
        #.reshape(ny_ccd,rebin,nx_ccd).sum(1)
        
        x_ccd_begin = int(floor(x_ccd_center-nx_ccd/2))
        y_ccd_begin = int(floor(y_ccd_center-ny_ccd/2))
        
        xx = slice(x_ccd_begin, (x_ccd_begin+nx_ccd))
        yy = slice(y_ccd_begin, (y_ccd_begin+ny_ccd))
        return xx,yy,ccd_pix_spot_values

    def _xypix_sincshift(self, ispec, wavelength):
        """
        Return xslice, yslice, pix for PSF at spectrum ispec, wavelength
        """
        print "in _xypix_sincshift"
        
        #- x,y of spot on CCD
        p, w = self._fiberpos[ispec], wavelength
        # xc = self._fx(p, w)
        # yc = self._fy(p, w)
        xc, yc = self.xy(ispec, wavelength)
        
        #- Ratio of CCD to Spot pixel sizes
        rpix = int(round(self.CcdPixelSize / self.SpotPixelSize))
        
        #- Calculate offset into CCD pixel
        xoffset = int(xc * rpix) % rpix
        yoffset = int(yc * rpix) % rpix

        #- Place high res spot into grid aligned with CCD pixels
        pix = self._fspot(p, w)
        ny, nx = pix.shape
        A = np.zeros(shape=(pix.shape[0]+rpix, pix.shape[1]+rpix))
        A[yoffset:yoffset+ny, xoffset:xoffset+nx] = pix
        ccdpix = rebin_image(A, rpix)
                
        #- Fractional high-res pixel offset
        #- This can be slow; is it really necessary?
        dxx = ((xc * rpix) % rpix - xoffset) / rpix
        dyy = ((yc * rpix) % rpix - yoffset) / rpix
        ccdpix = sincshift(ccdpix, dxx, dyy)
        
        #- sinc shift can cause negative ringing, so clip and re-normalize
        ccdpix = ccdpix.clip(0)
        ccdpix /= np.sum(ccdpix)

        #- Find where the [0,0] pixel goes on the CCD 
        #- Use floor() to get negative limits at edge of CCD correct
        xccd = int(np.floor(xc - ccdpix.shape[1]/2 + 1))
        yccd = int(np.floor(yc - ccdpix.shape[0]/2 + 1))
        
        xx = slice(xccd, xccd+ccdpix.shape[1])
        yy = slice(yccd, yccd+ccdpix.shape[0])
        
        return xx, yy, ccdpix

        
        
        
