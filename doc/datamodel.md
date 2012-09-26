Specter data model
==================

Input spectra
=============

### HDU 0 : blank ###

### HDU 1 : binary table of spectra ###

  - flux
  - calib
  - loglam or wavelength
  - sky
  - sky_calib

If `sky` exists, it is added to `flux`

If `calib` exists, `flux` and `sky` are in ergs/s/cm2/A and

    photons = (flux + sky) / calib
    
Otherwise, `flux` and `sky` are in units of photons as seen by CCD.

If `sky_calib` exists,

    photons = flux/calib + sky/sky_calib
    
This can be useful for simulating the throughput differences for 
pointing misalignments: the effective object throughput changes but not
the effective sky throughput.

i.e. depeneding upon whether `sky_calib`, `calib`, and `sky` exist,
the photons obeserved by the CCD could be:

    photons = flux/calib + sky/sky_calib
    photons = (flux + sky) / calib
    photons = flux / calib
    photons = flux + sky
    photons = flux


PSF Formats
===========

Base PSF
--------

HDU 0 : x[nspec, nwave]             EXTNAME="X"
HDU 1 : y[nspec, nwave]             EXTNAME="Y"
HDU 2 : wavelength[nspec, nwave]    EXTNAME="WAVELENGTH", or
        loglam[nspec, nwave]        EXTNAME="LOGLAM"
HDU 3+ : specific to each subformat

--> Additional extensions should be read by EXTNAME, not by number

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




























