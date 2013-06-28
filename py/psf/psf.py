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
from numpy.polynomial.legendre import Legendre
import scipy.optimize

from specter.util import gausspix
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
            
        #- Cache some ranges
        self._wmin = self._wavelength.min()
        self._wmax = self._wavelength.max()
            
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
        
        if wavelength < self._wavelength[ispec, 0] or \
           self._wavelength[ispec, -1] < wavelength:
            return slice(0,0), slice(0,0), N.zeros( (0,0) )
        
        xx, yy, ccdpix = self._xypix(ispec, wavelength)
        xlo, xhi = xx.start, xx.stop
        ylo, yhi = yy.start, yy.stop

        if xmax is None:
            xmax = self.npix_x
        if ymax is None:
            ymax = self.npix_y
        
        #- Check if completely off the edge in any direction
        if (xlo >= xmax) or (xhi <= xmin) or \
           (ylo >= ymax) or (yhi < ymin):
            return slice(0,0), slice(0,0), N.zeros( (0,0) )
            
        #- Check if partially off edge
        if xlo < xmin:
            ccdpix = ccdpix[:, -xlo:]
            xlo = xmin
        elif xhi > xmax:
            dx = xmax - xlo
            ccdpix = ccdpix[:, 0:dx]
            xhi = xmax

        if ylo < ymin:
            ccdpix = ccdpix[-ylo:, ]
            ylo = ymin
        elif yhi > ymax:
            dy = ymax - ylo
            ccdpix = ccdpix[0:dy, :]
            yhi = ymax
        
        xx = slice(xlo-xmin, xhi-xmin)
        yy = slice(ylo-ymin, yhi-ymin)
        
        return xx, yy, ccdpix

    def xyrange(self, spec_range, wave_range):
        """
        Return recommended range of pixels which cover these spectra/fluxes:
        (xmin, xmax, ymin, ymax)
        
        spec_range = indices specmin,specmax inclusive (not python style)
        wave_range = wavelength range wavemin,wavemax inclusive
        
        BUG: will fail if asking for a range where one of the spectra is
        completely off the CCD
        """
        specmin, specmax = spec_range
        wavemin, wavemax = wave_range

        #- Find the spectra with the smallest/largest y centroids
        ispec_ymin = specmin + N.argmin(self.y(None, wavemin)[specmin:specmax+1])
        ispec_ymax = specmin + N.argmax(self.y(None, wavemax)[specmin:specmax+1])
        ymin = self.xypix(ispec_ymin, wavemin)[1].start
        ymax = self.xypix(ispec_ymax, wavemax)[1].stop
        
        #- Now for wavelength where x = min(x),
        #- while staying on CCD and within wavelength range
        w = self.wavelength(specmin)
        ii = (wavemin < w) & (w < wavemax)
        ww = N.concatenate( ([wavemin,], w[ii], [wavemax,]) )
        x = self.x(specmin, ww)
        if min(x) < 0:
            xmin = 0.0
        else:
            wxmin = ww[N.argmin(x)]
            xmin = self.xypix(specmin, wxmin)[0].start
            
        #- and wavelength where x = max(x)
        w = self.wavelength(specmax)
        ii = (wavemin < w) & (w < wavemax)
        ww = N.concatenate( ([wavemin,], w[ii], [wavemax,]) )
        x = self.x(specmax, ww)
        if max(x) > self.npix_x:
            xmax = self.npix_x
        else:
            wxmax = ww[N.argmax(x)]
            xmax = self.xypix(specmax, wxmax)[0].stop
                                        
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
        
        ispec can be None or scalar
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
                
        if copy and isinstance(result, N.ndarray):
            return N.copy(result)
        else:
            return result

    def y(self, ispec=None, wavelength=None, copy=False):
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

        if copy and isinstance(result, N.ndarray):
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

    def xyw(self, ispec=None, copy=False):
        """
        Utility function to return x, y, and wavelength arrays for
        spectrum number ispec.  Set copy=True to get copies of the
        arrays instead of references to the originals.
        """
        x = self.x(ispec, copy=copy)
        y = self.y(ispec, copy=copy)
        w = self.wavelength(ispec, copy=copy)
        return x, y, w

    def loglam(self, ispec=None, y=None, copy=False):
        """
        Return log10(wavelength) of spectrum[ispec] evaluated at y.
        y can be None, scalar, or vector
        
        TODO: Directly interpolate on loglam grid?
        """
        return N.log10(self.wavelength(ispec, y, copy=copy))
    
    def wavelength(self, ispec=None, y=None, copy=False):
        """
        Return wavelength of spectrum[ispec] evaluated at y.
        
        ispec can be None or scalar; y can be None, scalar, or vector

        May return a view of the underlying array; do not modify unless
        specifying copy=True to get a copy of the data.
        """
        if ispec is None:
            if y is None:
                result = self._wavelength
            else:
                result = N.array([self.wavelength(i,y) for i in range(self.nspec)])
        else:
            if y is None:
                result = self._wavelength[ispec]
            else:
                result = N.interp(y, self._y[ispec], self._wavelength[ispec])
                
        if copy:
            return N.copy(result)
        else:
            return result
    
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

    def project(self, phot, wavelength, specmin=0, xr=None, yr=None, verbose=False):
        """
        Returns 2D image of spectra projected onto the CCD

        Required inputs:
            phot[nwave] or phot[nspec, nwave] as photons on CCD per bin
            wavelength[nwave] or wavelength[nspec, nwave] in Angstroms
                if wavelength is 1D and spectra is 2D, then wavelength[]
                applies to all phot[i]

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
                
            #- Only eval non-zero fluxes of wavelengths covered by this PSF
            wmin, wmax = self.wavelength(ispec, y=(0, self.npix_y))
            for j, w in enumerate(wspec):
                if phot[i,j] > 0.0 and wmin <= w and w <= wmax:
                    xx, yy, pix = self.xypix(ispec, w)
                    xx = slice(xx.start-xr[0], xx.stop-xr[0])
                    yy = slice(yy.start-yr[0], yy.stop-yr[0])
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
            spec_range = (ispecmin, ispecmax)
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
