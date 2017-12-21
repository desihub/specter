#!/usr/bin/env python
"""
SpotGridPSF - Linear interpolate hi-res sampled spots to model PSF

Stephen Bailey
Fall 2012
"""

from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import os
import numpy as np
from astropy.io import fits
from specter.psf import PSF
from specter.util import LinearInterp2D, rebin_image, sincshift
import scipy.interpolate
from numba import jit

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
        fx = fits.open(filename, memmap=False)
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
        return self._xypix_interp(ispec, wavelength)
        
    @jit(nopython=True,cache=True)
    def new_pixshift(dx,dy,w00,w10,w01,w11,pix_spot_values,rebin):
        """
        Return resampled_pix_spot_values 
        """
        ny_spot, nx_spot = pix_spot_values.shape
        #preallocate 
        resampled_pix_spot_values=np.zeros((ny_spot+rebin,nx_spot+rebin))
        for i in range(0,ny_spot):
            for j in range(0,nx_spot):  
                resampled_pix_spot_values[dy+i,dx+j]       += w00*pix_spot_values[i,j] 
                resampled_pix_spot_values[dy+1+i,dx+j]     += w10*pix_spot_values[i,j]
                resampled_pix_spot_values[dy+i,dx+1+j]     += w01*pix_spot_values[i,j]
                resampled_pix_spot_values[dy+1+i,dx+1+j]   += w11*pix_spot_values[i,j] 
                    
        return resampled_pix_spot_values
    
    def _xypix_interp(self, ispec, wavelength):
        """
        Return xslice, yslice, pix for PSF at spectrum ispec, wavelength
        """
        
       # # add itt tag
       # itt.resume()  
        
        #- Ratio of CCD to Spot pixel sizes
        rebin = int(self.CcdPixelSize / self.SpotPixelSize)
        
        p, w = self._fiberpos[ispec], wavelength
        pix_spot_values=self._fspot(p, w).astype(np.float64)
        nx_spot=pix_spot_values.shape[1]
        ny_spot=pix_spot_values.shape[0]
        nx_ccd=nx_spot//rebin+1 # add one bin because of resampling
        ny_ccd=ny_spot//rebin+1 # add one bin because of resampling
        
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
        resampled_pix_spot_values=new_pixshift(dx,dy,w00,w10,w01,w11,pix_spot_values,rebin)       
            
        # rebinning
        ccd_pix_spot_values=resampled_pix_spot_values.reshape(ny_spot+rebin,nx_ccd,rebin).sum(2).reshape(ny_ccd,rebin,nx_ccd).sum(1)
        # make sure it's positive
        ccd_pix_spot_values[ccd_pix_spot_values<0]=0.
        # normalize
        norm = np.sum(ccd_pix_spot_values)
        if norm > 0 :
            ccd_pix_spot_values /= norm

        x_ccd_begin = int(np.floor(xc))-nx_ccd//2+1  # begin of CCD coordinate stamp
        y_ccd_begin = int(np.floor(yc))-ny_ccd//2+1  # begin of CCD coordinate stamp
        xx = slice(x_ccd_begin, (x_ccd_begin+nx_ccd))
        yy = slice(y_ccd_begin, (y_ccd_begin+ny_ccd))
        

       # itt.detach()
        
        return xx,yy,ccd_pix_spot_values

        
        
        
    def _value(self, x, y, ispec, wavelength):

        """
        return PSF value (same shape as x and y), NOT integrated, for display of PSF.

        Arguments:
          x: x-coordinates baseline array
          y: y-coordinates baseline array (same shape as x)
          ispec: fiber
          wavelength: wavelength
        """

        
        p, w = self._fiberpos[ispec], wavelength
        pix_spot_values=self._fspot(p, w)
        nx_spot=pix_spot_values.shape[1]
        ny_spot=pix_spot_values.shape[0]
        x_spot=(np.arange(nx_spot)-nx_spot//2)
        y_spot=(np.arange(ny_spot)-ny_spot//2)
        spline=scipy.interpolate.RectBivariateSpline(x_spot,y_spot,pix_spot_values,kx=2, ky=2, s=0)
        
        xc, yc = self.xy(ispec, wavelength) # center of PSF in CCD coordinates
        ratio=self.CcdPixelSize/self.SpotPixelSize
        
        xr=x.ravel()
        yr=y.ravel()
        
        #img=spline((x-xc)*ratio,(y-yc)*ratio)
        #return img/np.sum(img)
        img=np.zeros(xr.shape)
        for i in range(xr.size) :
            img[i]=spline((xr[i]-xc)*ratio,(yr[i]-yc)*ratio)
        return img.reshape(x.shape)


        
        
        
