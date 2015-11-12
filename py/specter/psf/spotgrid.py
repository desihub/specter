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
        return self._xypix_interp(ispec, wavelength)
    
    def _xypix_interp(self, ispec, wavelength):
        """
        Return xslice, yslice, pix for PSF at spectrum ispec, wavelength
        """
        #- Ratio of CCD to Spot pixel sizes
        rebin = int(self.CcdPixelSize / self.SpotPixelSize)
        
        p, w = self._fiberpos[ispec], wavelength
        pix_spot_values=self._fspot(p, w)
        nx_spot=pix_spot_values.shape[1]
        ny_spot=pix_spot_values.shape[0]
        nx_ccd=nx_spot/rebin+1 # add one bin because of resampling
        ny_ccd=ny_spot/rebin+1 # add one bin because of resampling
        
        xc, yc = self.xy(ispec, wavelength) # center of PSF in CCD coordinates
                
        # fraction pixel offset requiring interpolation
        dx=xc*rebin-int(np.floor(xc*rebin)) # positive value between 0 and 1
        dy=yc*rebin-int(np.floor(yc*rebin)) # positive value between 0 and 1
        # weights for interpolation
        w00=(1-dy)*(1-dx)
        w10=dy*(1-dx)
        w01=(1-dy)*dx
        w11=dy*dx        
        # now the rest of the offset is an integer shift
        dx=int(np.floor(xc*rebin))-int(np.floor(xc))*rebin # positive integer between 0 and 14
        dy=int(np.floor(yc*rebin))-int(np.floor(yc))*rebin # positive integer between 0 and 14
        
        # resampled spot grid
        resampled_pix_spot_values=np.zeros((ny_spot+rebin,nx_spot+rebin))            
        resampled_pix_spot_values[dy:ny_spot+dy,dx:nx_spot+dx]         += w00*pix_spot_values
        resampled_pix_spot_values[dy+1:ny_spot+dy+1,dx:nx_spot+dx]     += w10*pix_spot_values
        resampled_pix_spot_values[dy:ny_spot+dy,dx+1:nx_spot+dx+1]     += w01*pix_spot_values
        resampled_pix_spot_values[dy+1:ny_spot+dy+1,dx+1:nx_spot+dx+1] += w11*pix_spot_values
            
        # rebinning
        ccd_pix_spot_values=resampled_pix_spot_values.reshape(ny_spot+rebin,nx_ccd,rebin).sum(2).reshape(ny_ccd,rebin,nx_ccd).sum(1)
        # make sure it's positive
        ccd_pix_spot_values[ccd_pix_spot_values<0]=0.
        # normalize
        n=np.sum(ccd_pix_spot_values)
        if n>0 :
            ccd_pix_spot_values /= n

        x_ccd_begin = int(np.floor(xc))-nx_ccd/2+1  # begin of CCD coordinate stamp
        y_ccd_begin = int(np.floor(yc))-ny_ccd/2+1  # begin of CCD coordinate stamp
        xx = slice(x_ccd_begin, (x_ccd_begin+nx_ccd))
        yy = slice(y_ccd_begin, (y_ccd_begin+ny_ccd))
        return xx,yy,ccd_pix_spot_values
        

    def _xypix_sincshift(self, ispec, wavelength):
        """
        Return xslice, yslice, pix for PSF at spectrum ispec, wavelength
        """
                
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

        
        
        
