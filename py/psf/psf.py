#!/usr/bin/env python

"""
Base class for 2D PSFs

Provides PSF base class which defines the interface for other code
using PSFs.  Subclasses implement specific models of the PSF and
override/extend the __init__ and xypix(ispec, wavelength) methods,
while allowing interchangeable use of different PSF models through
the interface defined in this base class.

Stephen Bailey, Fall 2012
"""

import numpy as N
### import scipy.sparse
import fitsio

class PSF(object):
    """
    Base class for 2D PSFs
    
    Subclasses need to extend __init__ to load format-specific items
    from the input fits file and implement xypix(ispec, wavelength)
    to return xslice, yslice, pixels[y,x] for the PSF evaluated at
    spectrum ispec at the given wavelength.  All interactions with PSF
    classes should be via the methods defined here, allowing
    interchangeable use of different PSF models
    """
    def __init__(self, filename):
        """
        Load PSF parameters from a file
        
        Loads x, y, wavelength information for spectral traces and fills:
            self.npix_x   #- number of columns in the target image
            self.npix_y   #- number of rows in the target image
            self.nspec    #- number of spectra (fibers)
            self.nwave    #- number of wavelength samples per spectrum

        Subclasses of this class define the xypix(ispec, wavelength) method
        to access the projection of this PSF into pixels.
        """
        
        #- Open fits file
        fx = fitsio.FITS(filename)
        hdr = fx[0].read_header()
        
        #- Load basic dimensions
        self.npix_x = hdr['NPIX_X']
        self.npix_y = hdr['NPIX_Y']
        self.nspec  = hdr['NAXIS2']
        self.nwave  = hdr['NAXIS1']
        
        #- Load x, y, loglam arrays
        self._x = fx[0].read()
        self._y = fx[1].read()
        if fx[2].get_extname() == 'WAVELENGTH':
            self._wavelength = fx[2].read()
            self._loglam = N.log10(self._wavelength)
        else:
            self._loglam = fx[2].read()
            self._wavelength = 10**self._loglam
        
    #-------------------------------------------------------------------------
    #- Evaluate the PSF into pixels
    
    def pix(self, ispec, wavelength):
        """
        Evaluate PSF for spectrum[ispec] at given wavelength
        
        returns 2D array pixels[iy,ix]
        
        also see xypix(ispec, wavelength)
        """
        return self.xypix(ispec, wavelength)[2]
        
    def xypix(self, ispec, wavelength):
        """
        Evaluate PSF for spectrum[ispec] at given wavelength
        
        returns xslice, yslice, pixels[iy,ix] such that
        image[yslice,xslice] += flux*pixels adds the contribution from
        spectrum ispec at that wavelength.
        
        Subclasses of PSF should implement this function for their
        specific model.
        """
        raise NotImplementedError

    def xyrange(self, spec_range, wave_range, dx=8, dy=8):
        """
        Return recommended range of pixels which cover these spectra/fluxes:
        (xmin, xmax, ymin, ymax)
        
        spec_range = indices specmin,specmax inclusive (not python style)
        wave_range = wavelength range wavemin,wavemax inclusive
        dx, dy = amount of extra overlap
        """
        specmin, specmax = spec_range
        wavemin, wavemax = wave_range

        #- Assume smallest x comes from specmin, and largest x from specmax
        xmin = N.min( self.x(specmin) )
        xmax = N.max( self.x(specmax) )

        #- Smallest/largest y could come from any spectrum
        yy = self._y[specmin:specmax+1]
        ww = self._wavelength[specmin:specmax+1]
        yy = yy[ N.where( (wavemin <= ww) & (ww <= wavemax) ) ]
        ymin = N.min(yy)
        ymax = N.max(yy)

        #- Make them integers
        xmin, xmax, ymin, ymax = map(int, map(round, (xmin, xmax, ymin, ymax)))
        
        #- Add borders, watching out for boundaries
        xmin = max(xmin-dx, 0)
        xmax = min(xmax+dx, self.npix_x)
        ymin = max(ymin-dy, 0)
        ymax = min(ymax+dy, self.npix_y)
        
        return (xmin, xmax, ymin, ymax)
    
    #-------------------------------------------------------------------------
    #- Shift PSF to a new x,y grid, e.g. to account for flexure
    
    def shift_xy(self, dx, dy):        
        """
        Shift the x,y trace locations of this PSF while preserving
        wavelength grid:  xnew = x + dx, ynew = y + dy
        """
        self._x += dx
        self._y += dy
    
    #-------------------------------------------------------------------------
    #- accessors for x, y, wavelength, loglam
        
    def x(self, ispec=None, wavelength=None, copy=False):
        """
        Return CCD X centroid of spectrum ispec at given wavelength(s).
        wavelength can be None, scalar, or a vector
        
        May return a view of the underlying array; do not modify unless
        specifying copy=True to get a copy of the data.
        """
        if ispec is None:
            if wavelength is None:
                result = self._x
            else:
                result = N.array([self.x(i, wavelength) for i in range(self.nspec)])
        else:
            if wavelength is None:
                result = self._x[ispec]
            else:
                result = N.interp(wavelength, self._wavelength[ispec], self._x[ispec])
                
        if copy:
            return N.copy(result)
        else:
            return result

    def y(self, ispec=None, wavelength=None, copy=True):
        """
        Return CCD Y centroid of spectrum ispec at given wavelength(s).
        wavelength can be None, scalar, or a vector

        May return a view of the underlying array; do not modify unless
        specifying copy=True to get a copy of the data.
        """
        if ispec is None:
            if wavelength is None:
                result = self._y
            else:
                result = N.array([self.y(i, wavelength) for i in range(self.nspec)])                
        else:
            if wavelength is None:
                result = self._y[ispec]
            else:
                result = N.interp(wavelength, self._wavelength[ispec], self._y[ispec])

        if copy:
            return N.copy(result)
        else:
            return result
            
    def xy(self, ispec=None, wavelength=None, copy=False):
        """
        Utility function to return self.x(...) and self.y(...) in one call
        """
        x = self.x(ispec, wavelength, copy=copy)
        y = self.y(ispec, wavelength, copy=copy)
        return x, y

    def loglam(self, ispec=None, y=None):
        """
        Return log10(wavelength) of spectrum[ispec] evaluated at y.
        y can be None, scalar, or vector
        """
        if ispec is None:
            return self._loglam
        elif y is None:
            return self._loglam[ispec]
        else:
            return N.interp(y, self._y[ispec], self._loglam[ispec])
    
    def wavelength(self, ispec=None, y=None):
        """
        Return wavelength of spectrum[ispec] evaluated at y.
        y can be None, scalar, or vector
        """
        if ispec is None:
            return self._wavelength
        elif y is None:
            return self._wavelength[ispec]
        else:
            return N.interp(y, self._y[ispec], self._wavelength[ispec])
    
    #-------------------------------------------------------------------------
    #- Project spectra onto CCD pixels
    def project(self, flux, wavelength, specmin=0, xr=None, yr=None):
        """
        Returns 2D image of spectra projected onto the CCD

        Required inputs:
            flux[nwave] or flux[nspec, nwave] as photons on CCD per bin
            wavelength[nwave] or wavelength[nspec, nwave] in Angstroms
                if wavelength is 1D and spectra is 2D, then wavelength[]
                applies to all flux[i]

        Optional inputs:
            specmin : starting spectrum number
            xr      : xrange xmin,xmax in CCD pixels
            yr      : yrange ymin,ymax in CCD pixels
        """
        #- x,y ranges and number of pixels
        if xr is None: xr = [0, self.npix_x]
        if yr is None: yr = [0, self.npix_y]
        nx = xr[1] - xr[0]    
        ny = yr[1] - yr[0]    

        #- For convenience, treat flux as a 2D vector
        flux = N.atleast_2d(flux)
        nspec, nw = flux.shape

        #- Create image to fill
        img = N.zeros( (ny, nx) )

        #- Loop over spectra and wavelengths
        for i, ispec in enumerate(range(specmin, specmin+nspec)):
            print ispec
            
            #- 1D wavelength for every spec, or 2D wavelength for 2D flux?
            if wavelength.ndim == 2:
                wspec = wavelength[i]
            else:
                wspec = wavelength
                
            #- Only eval non-zero fluxes of wavelengths covered by this PSF
            wpsf = self.wavelength(ispec)
            wmin, wmax = wpsf[0], wpsf[-1]
            for j, w in enumerate(wspec):
                if flux[i,j] > 0.0 and wmin <= w and w <= wmax:
                    xx, yy, pix = self.xypix(ispec, w)
                    xx = slice(xx.start-xr[0], xx.stop-xr[0])
                    yy = slice(yy.start-yr[0], yy.stop-yr[0])
                    img[yy, xx] += pix * flux[i,j]

        return img
    
    # #-------------------------------------------------------------------------
    # #- Access the projection matrix A
    # #- pix = A * flux
    #
    # #- NOTE: these are copied from bbspec PSF classes, used for extracting
    # #-       spectra, which isn't a part of specter (yet).
    # #-       This code is kept here for future reference.
    # 
    # def _reslice(self, xslice, yslice, data, xy_range):
    #     """
    #     shift xslice, yslice, and subsample data to match an xy_range,
    #     taking boundaries into account.
    #     """
    # 
    #     #- Check boundaries
    #     xmin, xmax, ymin, ymax = xy_range
    #     if xslice.start < xmin:
    #         dx = xmin - xslice.start
    #         xslice = slice(xmin, xslice.stop)
    #         data = data[:, dx:]
    #     if xslice.stop > xmax:
    #         dx = xslice.stop - xmax
    #         xslice = slice(xslice.start, xmax)
    #         data = data[:, :-dx]
    #     if yslice.start < ymin:
    #         dy = ymin - yslice.start
    #         yslice = slice(ymin, yslice.stop)
    #         data = data[dy:, :]
    #     if yslice.stop > ymax:
    #         dy = yslice.stop - ymax
    #         yslice = slice(yslice.start, ymax)
    #         data = data[:-dy, :]
    #         
    #     #- Shift slices    
    #     xslice = slice(xslice.start-xmin, xslice.stop-xmin)
    #     yslice = slice(yslice.start-ymin, yslice.stop-ymin)
    # 
    #     return xslice, yslice, data
    #     
    # 
    # def getAx(self, spec_range, flux_range, pix_range):
    #     """
    #     Returns sparse projection matrix from flux to pixels
    # 
    #     Inputs:
    #         spec_range = (ispecmin, ispecmax)
    #         flux_range = (ifluxmin, ifluxmax)
    #         pix_range  = (xmin, xmax, ymin, ymax)
    #     """
    # 
    #     #- Matrix dimensions
    #     specmin, specmax = spec_range
    #     fluxmin, fluxmax = flux_range
    #     xmin, xmax, ymin, ymax = pix_range        
    #     nspec = specmax - specmin
    #     nflux = fluxmax - fluxmin
    #     nx = xmax - xmin
    #     ny = ymax - ymin
    # 
    #     #- Generate A
    #     A = N.zeros( (ny*nx, nspec*nflux) )
    #     tmp = N.zeros((ny, nx))
    #     for ispec in range(specmin, specmax):
    #         for iflux in range(fluxmin, fluxmax):
    #             #- Get subimage and index slices
    #             xslice, yslice, pix = self.pix(ispec, iflux, xyrange=True)
    #             xslice, yslice, pix = self._reslice(xslice, yslice, pix, pix_range)
    #             
    #             #- If there is overlap with pix_range, put into sub-region of A
    #             if pix.shape[0]>0 and pix.shape[1]>0:              
    #                 tmp[yslice, xslice] = pix
    #                 ij = (ispec-specmin)*nflux + iflux-fluxmin
    #                 A[:, ij] = tmp.ravel()
    #                 tmp[yslice, xslice] = 0.0
    # 
    #     return scipy.sparse.csr_matrix(A)    
