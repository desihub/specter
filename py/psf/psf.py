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
from scipy.ndimage import center_of_mass
from numpy.polynomial.legendre import Legendre, legval, legfit
import scipy.optimize

from specter.util import gausspix, TraceSet
import fitsio

class PSF(object):
    """
    Base class for 2D PSFs
    
    Subclasses need to extend __init__ to load format-specific items
    from the input fits file and implement _xypix(ispec, wavelength)
    to return xslice, yslice, pixels[y,x] for the PSF evaluated at
    spectrum ispec at the given wavelength.  All interactions with PSF
    classes should be via the methods defined here, allowing
    interchangeable use of different PSF models.
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
        
        #- Load basic dimensions
        hdr = fitsio.read_header(filename)
        self.npix_x = hdr['NPIX_X']
        self.npix_y = hdr['NPIX_Y']
        self.nspec  = hdr['NSPEC']
        
        #- Load x, y legendre coefficient tracesets
        xc, hdr = fitsio.read(filename, 'XCOEFF', header=True)
        self._x = TraceSet(xc, domain=(hdr['WAVEMIN'], hdr['WAVEMAX']))
        yc, hdr = fitsio.read(filename, 'YCOEFF', header=True)
        self._y = TraceSet(yc, domain=(hdr['WAVEMIN'], hdr['WAVEMAX']))
        
        #- Create inverse y -> wavelength mapping
        self._w = self._y.invert()
        self._wmin = N.min(self.wavelength(None, 0))
        self._wmax = N.max(self.wavelength(None, self.npix_y-1))
                
        #- Filled only if needed
        self._xsigma = None
        self._ysigma = None

    #- Utility function to fit spot sigma vs. wavelength
    def _fit_spot_sigma(self, ispec, axis=0, npoly=5):
        """
        Fit the cross-sectional Gaussian sigma of PSF spots vs. wavelength.
        Return callable Legendre object.
        
        Inputs:
            ispec : spectrum number
            axis  : 0 or 'x' for cross dispersion sigma;
                    1 or 'y' or 'w' for wavelength dispersion
            npoly : order of Legendre poly to fit to sigma vs. wavelength
            
        Returns:
            legfit such that legfit(w) returns fit at wavelengths w
        """
        
        if type(axis) is not int:
            if axis in ('x', 'X'):
                axis = 0
            elif axis in ('y', 'Y', 'w', 'W'):
                axis = 1
            else:
                raise ValueError("Unknown axis type "+str(axis))
                
        if axis not in (0,1):
            raise ValueError("axis must be 0, 'x', 1, 'y', or 'w'")
            
        yy = N.linspace(10, self.npix_y-10, 20)
        ww = self.wavelength(ispec, y=yy)
        xsig = list()  #- sigma vs. wavelength array to fill
        for w in ww:
            xspot = self.pix(ispec, w).sum(axis=axis)
            xspot /= N.sum(xspot)       #- normalize for edge cases
            xx = N.arange(len(xspot))
            mean, sigma = scipy.optimize.curve_fit(gausspix, xx, xspot)[0]
            xsig.append(sigma)
        
        #- Fit Legendre polynomial and return coefficients
        legfit = Legendre.fit(ww, xsig, npoly, domain=(self._wmin, self._wmax))
                
        return legfit

    #-------------------------------------------------------------------------
    #- Cross dispersion width for row-by-row extractions
    def xsigma(self, ispec, wavelength):
        """
        Return Gaussian sigma of PSF spot in cross-dispersion direction
        in CCD pixel units.
        
        ispec : spectrum index
        wavelength : scalar or vector wavelength(s) to evaluate spot sigmas
        
        The first time this is called for a spectrum, the PSF is sampled
        at 20 wavelengths and the variation is fit with a 5th order
        Legendre polynomial and the coefficients are cached.
        The actual value (and subsequent calls) use these cached
        Legendre fits to interpolate the sigma value.  If this is not
        fast enough and/or accurate enough, PSF subtypes may override
        this function to provide a more accurate xsigma measurement.
        """

        #- First call for any spectrum: setup array to cache coefficients
        if self._xsigma is None:
            self._xsigma = [None,] * self.nspec
            
        #- First call for this spectrum: calculate coefficients & cache
        if self._xsigma[ispec] is None:
            self._xsigma[ispec] = self._fit_spot_sigma(ispec, axis=0, npoly=5)

        #- Use cached Legendre fit to interpolate xsigma at wavelength(s)
        return self._xsigma[ispec](wavelength)

    #-------------------------------------------------------------------------
    #- Cross dispersion width for row-by-row extractions
    def ysigma(self, ispec, wavelength):
        """
        Return Gaussian sigma of PSF spot in wavelength-dispersion direction
        in units of pixels.
        
        Also see wdisp(...) which returns sigmas in units of Angstroms.
        
        ispec : spectrum index
        wavelength : scalar or vector wavelength(s) to evaluate spot sigmas
        
        See notes in xsigma(...) about caching of Legendre fit coefficients.
        """

        #- First call for any spectrum: setup array to cache coefficients
        if self._ysigma is None:
            self._ysigma = [None,] * self.nspec
            
        #- First call for this spectrum: calculate coefficients & cache
        if self._ysigma[ispec] is None:
            self._ysigma[ispec] = self._fit_spot_sigma(ispec, axis=1, npoly=5)

        #- Use cached Legendre fit to interpolate xsigma at wavelength(s)
        return self._ysigma[ispec](wavelength)

    #-------------------------------------------------------------------------
    #- Cross dispersion width for row-by-row extractions
    def wdisp(self, ispec, wavelength):
        """
        Return Gaussian sigma of PSF spot in wavelength-dispersion direction
        in units of Angstroms.
        
        Also see ysigma(...) which returns sigmas in units of pixels.
        
        ispec : spectrum index
        wavelength : scalar or vector wavelength(s) to evaluate spot sigmas
        
        See notes in xsigma(...) about caching of Legendre fit coefficients.
        """
        sigma_pix = self.ysigma(ispec, wavelength)
        return self.angstroms_per_pixel(ispec, wavelength) * sigma_pix        

    #-------------------------------------------------------------------------
    #- Evaluate the PSF into pixels
    
    def pix(self, ispec, wavelength):
        """
        Evaluate PSF for spectrum[ispec] at given wavelength
        
        returns 2D array pixels[iy,ix]
        
        also see xypix(ispec, wavelength)
        """
        return self.xypix(ispec, wavelength)[2]

    def _xypix(self, ispec, wavelength):
        """
        Subclasses of PSF should implement this to return
        xslice, yslice, pixels[iy,ix] for their particular
        models.  Don't worry about edge effects -- PSF.xypix
        will take care of that.
        """
        raise NotImplementedError
        
    def xypix(self, ispec, wavelength, xmin=0, xmax=None, ymin=0, ymax=None):
        """
        Evaluate PSF for spectrum[ispec] at given wavelength
        
        returns xslice, yslice, pixels[iy,ix] such that
        image[yslice,xslice] += photons*pixels adds the contribution from
        spectrum ispec at that wavelength.
        
        if xmin or ymin are set, the slices are relative to those
        minima (useful for simulating subimages)
        """
        
        if xmax is None:
            xmax = self.npix_x
        if ymax is None:
            ymax = self.npix_y

        if wavelength < self.wavelength(ispec, -0.5):
            return slice(0,0), slice(0,0), N.zeros((0,0))
        elif wavelength > self.wavelength(ispec, self.npix_y-0.5):
            return slice(0,0), slice(ymax, ymax), N.zeros((0,0))
        
        xx, yy, ccdpix = self._xypix(ispec, wavelength)
        xlo, xhi = xx.start, xx.stop
        ylo, yhi = yy.start, yy.stop

        #- Check if completely off the edge in any direction
        if (ylo >= ymax):
            return slice(0,0), slice(ymax,ymax), N.zeros( (0,0) )            
        elif (yhi < ymin):
            return slice(0,0), slice(ymin,ymin), N.zeros( (0,0) )
        elif (xlo >= xmax):
            return slice(xmax, xmax), slice(0,0), N.zeros( (0,0) )
        elif (xhi <= xmin):
            return slice(xmin, xmin), slice(0,0), N.zeros( (0,0) )
            
        #- Check if partially off edge
        if xlo < xmin:
            ccdpix = ccdpix[:, -(xhi-xmin):]
            xlo = xmin
        elif xhi > xmax:
            ccdpix = ccdpix[:, 0:(xmax-xlo)]
            xhi = xmax

        if ylo < ymin:
            ccdpix = ccdpix[-(yhi-ymin):, ]
            ylo = ymin
        elif yhi > ymax:
            ccdpix = ccdpix[0:(ymax-ylo), :]
            yhi = ymax
        
        xx = slice(xlo-xmin, xhi-xmin)
        yy = slice(ylo-ymin, yhi-ymin)
        
        #- Check if we are off the edge
        if (xx.stop-xx.start == 0) or (yy.stop-yy.start == 0):
            ccdpix = N.zeros( (0,0) )
                
        return xx, yy, ccdpix

    def xyrange(self, spec_range, wavelengths):
        """
        Return recommended range of pixels which cover these spectra/fluxes:
        (xmin, xmax, ymin, ymax)
        
        spec_range = indices specmin,specmax (python style indexing),
                     or scalar for single spectrum index
        wavelengths = wavelength range wavemin,wavemax inclusive
                     or sorted array of wavelengths
        
        BUG: will fail if asking for a range where one of the spectra is
        completely off the CCD
        """
        if isinstance(spec_range, int):
            specmin, specmax = spec_range, spec_range+1
        else:
            specmin, specmax = spec_range

        if isinstance(wavelengths, (int, float)):
            wavemin = wavemax = wavelengths
        else:
            wavemin, wavemax = wavelengths[0], wavelengths[-1]

        #- Find the spectra with the smallest/largest y centroids
        ispec_ymin = specmin + N.argmin(self.y(None, wavemin)[specmin:specmax+1])
        ispec_ymax = specmin + N.argmax(self.y(None, wavemax)[specmin:specmax+1])
        ymin = self.xypix(ispec_ymin, wavemin)[1].start
        ymax = self.xypix(ispec_ymax, wavemax)[1].stop
        
        #- Now for wavelength where x = min(x),
        #- while staying on CCD and within wavelength range
        w = self.wavelength(specmin)
        if w[0] < wavemin:
            w = w[wavemin <= w]
        if wavemax < w[-1]:
            w = w[w <= wavemax]
        
        #- Add in wavemin and wavemax since w isn't perfect resolution
        w = N.concatenate( (w, (wavemin, wavemax) ) )
        
        x = self.x(specmin, w)
        if min(x) < 0:
            xmin = 0.0
        else:
            wxmin = w[N.argmin(x)]  #- wavelength at x minimum
            xmin = self.xypix(specmin, wxmin)[0].start
            
        #- and wavelength where x = max(x)
        w = self.wavelength(specmax-1)
        if w[0] < wavemin:
            w = w[wavemin <= w]
        if wavemax < w[-1]:
            w = w[w <= wavemax]
        
        #- Add in wavemin and wavemax since w isn't perfect resolution
        w = N.concatenate( (w, (wavemin, wavemax) ) )

        x = self.x(specmax-1, w)
        if max(x) > self.npix_x:
            xmax = self.npix_x
        else:
            wxmax = w[N.argmax(x)]
            xmax = self.xypix(specmax-1, wxmax)[0].stop
                                        
        return (xmin, xmax, ymin, ymax)
    
    #-------------------------------------------------------------------------
    #- Shift PSF to a new x,y grid, e.g. to account for flexure
    
    def shift_xy(self, dx, dy):
        """
        Shift the x,y trace locations of this PSF while preserving
        wavelength grid:  xnew = x + dx, ynew = y + dy
        """
        raise NotImplementedError
    
    #-------------------------------------------------------------------------
    #- accessors for x, y, wavelength
                
    def x(self, ispec=None, wavelength=None):
        """
        Return CCD X centroid of spectrum ispec at given wavelength(s).
        
        ispec can be None, scalar, or vector
        wavelength can be None, scalar or a vector
        
        ispec   wavelength  returns
        +-------+-----------+------
        None    None        array[nspec, npix_y]
        None    scalar      vector[nspec]
        None    vector      array[nspec, nwave]
        scalar  None        array[npix_y]
        scalar  scalar      scalar
        scalar  vector      vector[nwave]
        vector  None        array[nspec, npix_y]
        vector  scalar      vector[nspec]
        vector  vector      array[nspec, nwave]
        """
        
        if wavelength is None:
            #- ispec=None -> ispec=every spectrum
            if ispec is None:
                ispec = N.arange(self.nspec)
            
            #- ispec is an array; sample at every row
            if isinstance(ispec, (N.ndarray, list, tuple)):
                x = list()
                for i in ispec:
                    w = self.wavelength(i)
                    x.append(self._x.eval(i, w))
                return N.array(x)
            else:  #- scalar ispec, make wavelength an array
                wavelength = self.wavelength(ispec)
            
        return self._x.eval(ispec, wavelength)
            
    def y(self, ispec=None, wavelength=None):
        """
        Return CCD Y centroid of spectrum ispec at given wavelength(s).
        
        ispec can be None, scalar, or vector
        wavelength can be scalar or a vector (but not None)
        
        ispec   wavelength  returns
        +-------+-----------+------
        None    scalar      vector[nspec]
        None    vector      array[nspec,nwave]
        scalar  scalar      scalar
        scalar  vector      vector[nwave]
        vector  scalar      vector[nspec]
        vector  vector      array[nspec, nwave]
        """
        if wavelength is None:
            raise ValueError, "PSF.y requires wavelength scalar or vector"
            
        if ispec is None:
            ispec = N.arange(self.nspec)
            
        return self._y.eval(ispec, wavelength)
        
        if ispec is None:
            if wavelength is None:
                return N.tile(N.arange(self.npix_y), self.nspec).reshape(self.nspec, self.npix_y)
            else:
                ispec = N.arange(self.nspec)
        
        if wavelength is None:
            wavelength = self.wavelength(ispec)

        return self._y.eval(ispec, wavelength)
            
    def xy(self, ispec=None, wavelength=None):
        """
        Utility function to return self.x(...) and self.y(...) in one call
        """
        x = self.x(ispec, wavelength)
        y = self.y(ispec, wavelength)
        return x, y

    def wavelength(self, ispec=None, y=None):
        """
        Return wavelength of spectrum[ispec] evaluated at y.
        
        ispec can be None, scalar, or vector
        y can be None, scalar, or vector

        May return a view of the underlying array; do not modify unless
        specifying copy=True to get a copy of the data.
        """
        if y is None:
            y = N.arange(0, self.npix_y)
        
        if ispec is None:
            ispec = N.arange(self.nspec)
            
        return self._w.eval(ispec, y)
    
    def angstroms_per_pixel(self, ispec, wavelength):
        """
        Return CCD pixel width in Angstroms for spectrum ispec at given
        wavlength(s).  Wavelength may be scalar or array.
        """
        ww = self.wavelength(ispec, y=N.arange(self.npix_y))
        dw = N.gradient( ww )
        return N.interp(wavelength, ww, dw)
    
    #-------------------------------------------------------------------------
    #- Project spectra onto CCD pixels
    # def project_subimage(self, phot, wavelength, specmin, verbose=False):
    #     """
    #     Project photons onto CCD.  Returns subimage, (xmin,xmax,ymin,ymax).
    #     See PSF.project() for full parameter descriptions.
    #     """
    #     #- NOTES:
    #     #- Tightly coupled to self.project
    #     #- Should this return slices instead of xyrange, similar to
    #     #-     PSF.xypix?
    #     #- Maybe even rename to xyproject() ?
    #     
    #     nspec = phot.shape[0] if phot.ndim == 2 else self.nspec
    #     specmax = min(specmin+nspec, nspec)
    #     specrange = (specmin, specmax)
    #     waverange = (N.min(wavelength), N.max(wavelegth))
    #     xmin, xmax, ymin, ymax = xyrange = self.xyrange(specrange, waverange)
    #     image = self.project(phot, wavelength, specmin=specmin, \
    #         xr=(xmin,xmax), yr=(ymin, ymax), verbose=verbose)
    #         
    #     return image, xyrange

    def project(self, phot, wavelength, specmin=0, xyrange=None, verbose=False):
        """
        Returns 2D image of spectra projected onto the CCD

        Required inputs:
            phot[nwave] or phot[nspec, nwave] as photons on CCD per bin
            wavelength[nwave] or wavelength[nspec, nwave] in Angstroms
                if wavelength is 1D and spectra is 2D, then wavelength[]
                applies to all phot[i]

        Optional inputs:
            specmin : starting spectrum number
            xyrange : (xmin, xmax, ymin, ymax) range of CCD pixels
        """
        #- x,y ranges and number of pixels
        if xyrange is None:
            xmin, xmax = (0, self.npix_x)
            ymin, ymax = (0, self.npix_y)
            xyrange = (xmin, xmax, ymin, ymax)
        else:
            xmin, xmax, ymin, ymax = xyrange
            
        nx = xmax - xmin
        ny = ymax - ymin

        #- For convenience, treat phot as a 2D vector
        phot = N.atleast_2d(phot)
        nspec, nw = phot.shape

        #- Create image to fill
        img = N.zeros( (ny, nx) )

        #- Loop over spectra and wavelengths
        for i, ispec in enumerate(range(specmin, specmin+nspec)):
            if verbose:
                print ispec
            
            #- 1D wavelength for every spec, or 2D wavelength for 2D phot?
            if wavelength.ndim == 2:
                wspec = wavelength[i]
            else:
                wspec = wavelength
                
            #- Evaluate positive photons within wavelength range
            wmin, wmax = self.wavelength(ispec, y=(0, self.npix_y))
            for j, w in enumerate(wspec):
                if phot[i,j] > 0.0 and (wmin <= w <= wmax):
                    xx, yy, pix = self.xypix(ispec, w, \
                        xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
                    img[yy, xx] += pix * phot[i,j]

        return img
    
    #- Convenience functions
    
    @property
    def wmin(self):
        return self._wmin

    @property
    def wmax(self):
        return self._wmax
    
    def projection_matrix(self, spec_range, wavelengths, xyrange):
        """
        Returns sparse projection matrix from flux to pixels
    
        Inputs:
            spec_range = (ispecmin, ispecmax) or scalar ispec
            wavelengths = array_like wavelengths
            xyrange  = (xmin, xmax, ymin, ymax)
            
        Usage:
            xyrange = xmin, xmax, ymin, ymax
            A = psf.projection_matrix(spec_range, wavelengths, xyrange)
            nx = xmax-xmin
            ny = ymax-ymin
            img = A.dot(phot.ravel()).reshape((ny,nx))
        """
    
        #- Matrix dimensions
        if isinstance(spec_range, int):
            specmin, specmax = spec_range, spec_range+1
        else:
            specmin, specmax = spec_range
            
        xmin, xmax, ymin, ymax = xyrange        
        nspec = specmax - specmin
        nflux = len(wavelengths)
        nx = xmax - xmin
        ny = ymax - ymin
    
        #- Generate A
        A = N.zeros( (ny*nx, nspec*nflux) )
        tmp = N.zeros((ny, nx))
        for ispec in range(specmin, specmax):
            for iflux, w in enumerate(wavelengths):
                #- Get subimage and index slices
                xslice, yslice, pix = self.xypix(ispec, w, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
                
                #- If there is overlap with pix_range, put into sub-region of A
                if pix.shape[0]>0 and pix.shape[1]>0:         
                    tmp[yslice, xslice] = pix
                    ij = (ispec-specmin)*nflux + iflux
                    A[:, ij] = tmp.ravel()
                    tmp[yslice, xslice] = 0.0
    
        return scipy.sparse.csr_matrix(A)    

