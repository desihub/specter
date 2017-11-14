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
import time

class SpotGridPSF(PSF):
    """
    Model PSF with a linear interpolation of high resolution sampled spots
    """
    xypix_interp_elapsed_t=0
    
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
    
    def _xypix_interp(self, ispec, wavelength):
        """
        Return xslice, yslice, pix for PSF at spectrum ispec, wavelength
        """
        #add timer for whole function ---------------------------------------------
        xypix_interp_t0=time.time()
        
        #- Ratio of CCD to Spot pixel sizes

        rebin = int(self.CcdPixelSize / self.SpotPixelSize)
        
        #add timer for fiberpos ----------------------
        fiberpos_t0=time.time()
        p, w = self._fiberpos[ispec], wavelength
        fiberpos_t1=time.time()
        #done timing fiberpos -----------------------
        fiberpos_elapsed_t=fiberpos_t1-fiberpos_t0
        
        #add timer for fspot -------------------------
        #this is the slower one (20 percent), fspot calls LinearInterp2D, which is probably the slow part
        fspot_t0=time.time()
        pix_spot_values=self._fspot(p, w)
        fspot_t1=time.time()
        #done timing fspot ---------------------------
        fspot_elapsed_t=fspot_t1-fspot_t0
        
        nx_spot=pix_spot_values.shape[1]
        ny_spot=pix_spot_values.shape[0]
        nx_ccd=nx_spot//rebin+1 # add one bin because of resampling
        ny_ccd=ny_spot//rebin+1 # add one bin because of resampling
        
        #add timer for xy----------------------------
        xy_t0=time.time()
        xc, yc = self.xy(ispec, wavelength) # center of PSF in CCD coordinates
        xy_t1=time.time()
        #done timing for xy ------------------------
        xy_elapsed_t=xy_t1-xy_t0
                
        #timer for pixel offset --------------------
        offset_t0=time.time()
                
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
        
        offset_t1=time.time()
        offset_elapsed_t=offset_t1-offset_t0
        #done timing offset -----------------------
        #print("offset elapsed time is %s s" %(offset_elapsed_t))
        
        # resampled spot grid     
        #start timer for zeros creation ------------
        zeros_t0=time.time()        
        resampled_pix_spot_values=np.zeros((ny_spot+rebin,nx_spot+rebin))     
       
        zeros_t1=time.time()
        #done timing zeros -------------------
        zeros_elapsed_t=zeros_t1-zeros_t0

        #start timer for resampling grid -----------------------
        resample_t0=time.time()    

        resampled_pix_spot_values[dy:ny_spot+dy,dx:nx_spot+dx]         += w00*pix_spot_values
        resampled_pix_spot_values[dy+1:ny_spot+dy+1,dx:nx_spot+dx]     += w10*pix_spot_values
        resampled_pix_spot_values[dy:ny_spot+dy,dx+1:nx_spot+dx+1]     += w01*pix_spot_values
        resampled_pix_spot_values[dy+1:ny_spot+dy+1,dx+1:nx_spot+dx+1] += w11*pix_spot_values
        
        resample_t1=time.time()
        #done timing resample ------------------------------
        resample_elapsed_t=resample_t1-resample_t0
        #print("resample elapsed time is %s s" %(resample_elapsed_t))
            
        #start timing ccd_rebin -------------------------------
        ccd_rebin_t0=time.time()
        # rebinning
        ccd_pix_spot_values=resampled_pix_spot_values.reshape(ny_spot+rebin,nx_ccd,rebin).sum(2).reshape(ny_ccd,rebin,nx_ccd).sum(1)
        # make sure it's positive
        ccd_pix_spot_values[ccd_pix_spot_values<0]=0.
        # normalize
        norm = np.sum(ccd_pix_spot_values)
        if norm > 0 :
            ccd_pix_spot_values /= norm
            
        ccd_rebin_t1=time.time()
        #done timing ccd_rebin ----------------------------------
        ccd_rebin_elapsed_t=ccd_rebin_t1-ccd_rebin_t0
        #print("ccd_rebin elapsed time is % s" %(ccd_rebin_elapsed_t))
        
        #start timing ccd_slice ------------------------------
        ccd_slice_t0=time.time()

        x_ccd_begin = int(np.floor(xc))-nx_ccd//2+1  # begin of CCD coordinate stamp
        y_ccd_begin = int(np.floor(yc))-ny_ccd//2+1  # begin of CCD coordinate stamp
        xx = slice(x_ccd_begin, (x_ccd_begin+nx_ccd))
        yy = slice(y_ccd_begin, (y_ccd_begin+ny_ccd))
        
        ccd_slice_t1=time.time()
        #done timing ccd_slice -------------------------------
        ccd_slice_elapsed_t=ccd_slice_t1-ccd_slice_t0
        #print("ccd_slice elapsed time is % s" %(ccd_slice_elapsed_t))
        

        #add final timer
        xypix_interp_t1=time.time()
        #add previously stored value, initially zero
        xypix_interp_elapsed_t=xypix_interp_elapsed_t + xypix_interp_t1-xypix_interp_t0
        
        #done timing -------------------------------------------------------
        #print("_xypix_interp elapsed time is %s s" %(xypix_elapsed_t))
        
        #now compute fraction of time each part of xypix_interp each time block takes
        fiberpos_frac=fiberpos_elapsed_t/xypix_interp_elapsed_t
        fspot_frac=fspot_elapsed_t/xypix_interp_elapsed_t
        xy_frac=xy_elapsed_t/xypix_interp_elapsed_t
        offset_frac=offset_elapsed_t/xypix_interp_elapsed_t
        zeros_frac=zeros_elapsed_t/xypix_interp_elapsed_t
        resample_frac=resample_elapsed_t/xypix_interp_elapsed_t
        ccd_rebin_frac=ccd_rebin_elapsed_t/xypix_interp_elapsed_t
        ccd_slice_frac=ccd_slice_elapsed_t/xypix_interp_elapsed_t
        #for a sanity check, check total fraction tracked
        xypix_interp_frac=fiberpos_frac + fspot_frac + xy_frac + offset_frac + zeros_frac + resample_frac + ccd_rebin_frac + ccd_slice_frac
        
        #print("xypix_interp fiberpos fraction used is %s" %(fiberpos_frac))
        #print("xypix_interp fspot fraction used is %s" %(fspot_frac))
        #print("xypix_interp xy fraction used is %s" %(xy_frac))
        #print("xypix_interp offset fraction used is %s" %(offset_frac))
        #print("xypix_interp zeros fraction used is %s" %(zeros_frac))
        #print("xypix_interp resample fraction used is %s" %(resample_frac))
        #print("xypix_interp ccd_rebin fraction used is %s" %(ccd_rebin_frac))
        #print("xypix_interp ccd_slice fraction used is %s" %(ccd_slice_frac))
        #print("total xypix_interp tracked is %s"%(xypix_interp_frac))
        
        print("runtime for xypix_interp is %s s" %(xypix_interp_elapsed_t))
        
        #return the elapsed time in the function so we can continue aggregating statistics
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


        
        
        
