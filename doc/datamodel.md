Specter data model
==================

Input spectra
=============

### HDU 0 : blank ###

### HDU 1 : binary table of spectra ###

Required:
  - flux[nwave] or flux[nspec, nwave]
  - loglam or wavelength (1D[nwave] or 2D[nspec, nwave])

Optional:
  - sky[nwave] or sky[nspec, nwave]

If `sky` exists, it is added to `flux`.  Other user-specific columns may
exist and will be ignored by specter.

If units aren't specified, flux will be treated as photons/bin as seen
by the CCD.  If TUNITn keyword for the flux column starts with "erg" then
it will be treated as ergs/sec/cm^2/A and converted to photons and throughput
applied.  In both cases, the flux will be treated as a delta function at
that wavelength; it will not be integrated or interpolated between points.
Thus the wavelength sampling should be finer than the PSF resolution.


PSF Formats
===========

Base PSF
--------

All PSF files have the same format for HDUs 0-2; additional extensions
are specific to each subtype and contain whatever information in whatever
format is needed to express the PSF in that parameterization.

HDU 0 : x[nspec, nwave]             EXTNAME="X"
HDU 1 : y[nspec, nwave]             EXTNAME="Y"
HDU 2 : wavelength[nspec, nwave]    EXTNAME="WAVELENGTH", or
        loglam[nspec, nwave]        EXTNAME="LOGLAM"
HDU 3+ : specific to each subformat

--> Additional extensions should be read by EXTNAME, not by number; the
    order of these extensions should be treated as arbitrary

HDU 0 keywords:
  - NPIX_X, NPIX_Y : CCD dimensions in pixels

Optional: extensions and keywords to convert ergs/s/cm^2/A -> photons

HDU EXTNAME="THROUGHPUT", binary table with columns:
  - wavelength[nwave] or loglam[nwave]
  - throughput[nwave]

HDU 0 keywords
  - EXPTIME: Standard exposure time in seconds
  - EFFAREA: Effective area of mirror, including obscuration effects, in cm^2
  
photons/A = (ergs/s/cm^2/A) * EXPTIME * EFFAREA * wavelength / (h*c)

These these keywords aren't available, the PSF can still be used to project
photons onto a CCD, but not flux in ergs/s/cm^2/A .

Spot Grid PSF
-------------
PSFTYPE = "SPOTGRID"

This PSF type provides an x,y grid of spots to be interpolated.
The first interpolation dimension must monotonically increase with
spectrum number, e.g. the position of a fiber along a slit head.
The second interpolation dimension is wavelength.

HDU 0-2 : Same as Base PSF: X, Y, wavelength or loglam of traces

HDU SPOTS : spot[nspotx, nspoty, iy, ix]    #- 2D PSF spots
    NAXIS1 = number of spot samples in the spectrum number dimension
    NAXIS2 = number of spot samples in the wavelength direction
    NAXIS3 = size of spot in the CCD y dimension
    NAXIS4 = size of spot in the CCD x direction
    
HDU XCCD : xccd[nspotx, nspoty]     #- CCD X pixel location of spot pixel[0,0]
HDU YCCD : yccd[nspotx, nspoty]     #- CCD Y pixel location of spot pixel[0,0]
HDU FIBERPOS : fiberpos[nspec]      #- Slit position of each fiber
HDU SPOTPOS  : spotpos[nspotx]      #- Slit positions where spots are sampled
HDU SPOTWAVE : spotwave[nspoty]     #- Wavelengths where spots are sampled

spot[i,j] is a 2D PSF spot sampled at slit position spotpos[i] and
wavelength spotwave[j].  Its [0,0] pixel is located on the CCD at
xccd[i,j], yccd[i,j].




























