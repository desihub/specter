Specter data model
==================

Input Spectra
=============

### HDU 0 : blank ###

### HDU 1 : binary table of spectra ###

Required:
  - flux[nwave] or flux[nspec, nwave]
  - header keywords CRVAL1 and CDELT1, or
    `wave`, `wavelength` or `loglam` column with dimensions matching `flux`

Optional:
  - objtype[nspec] column, or header keyword OBJTYPE.  Default 'STAR'.
    Values = CALIB, SKY, or something else.  Anything other
    than CALIB or SKY will be treated as an astronomical source.
    See the Throughput section for a description of how the throughput
    model is applied differently for various object types.

Other columns and HDUs may be present and will be ignored.

### Wavelength Grid ###

The spectral wavelengths can be specified with either header keywords
CRVAL1 and CDELT1 and optionally DC-FLAG, or with an arbitrary wavelength
grid in a `wavelength` or `loglam` column whose dimensions match `flux`.
DC-FLAG = 0/1 for linear/log10 wavelength spacing; default=0.

e.g. to specify wavelengths [3600, 3601, 3602, ...]:

    CRVAL1  =                 3600 / Reference wavelength in A
    CDELT1  =                    1 / delta wavelength
    DC-FLAG =                    0 / Linearly spaced wavelengths

or to specify log10(wavelength) spacing [3.5500, 3.5501, 3.5502, ...]

    CRVAL1  =               3.5500 / Reference log10(Angstroms)
    CDELT1  =               0.0001 / delta log10(Angstroms)
    DC-FLAG =                    0 / Log10 spaced wavelengths

The FITS standard also requires CRPIX1 and CTYPE1 but these are ignored.

### Flux Units ###

If the `flux` column has a TUNITnn keyword, that will be used, otherwise
keyword FLUXUNIT is an alternative.  Options:
  * Treated as function values to be multipled by bin width:
    - erg/s/cm^2/A  (default)
    - erg/s/cm^2/A/arcsec^2
    - photon/A
    - photon/A/arcsec^2
  * Treated as delta functions at each given wavelength:
    - photons
    - erg/s/cm^2

For example, an astromical object is typically in units "erg/s/cm^2/A"
and will be converted to photons using all throughput terms of the
throughput model.  A sky spectrum may be in "erg/s/cm^2/A/arcsec^2" and
will be multiplied by the area of the fiber instead of having a
fiber input geometric loss applied.


Throughput
==========

### HDU EXTNAME=THROUGHPUT ###

binary table with columns:
  - wavelength[nwave]   Wavelength in Angstroms,
      or loglam[nwave]  for log10(Angstroms)
  - extinction[nwave]   Atmospheric extinction in mags/airmass;
                        Atm throughput = 10^(0.4 * extinction * airmass)
  - throughput[nwave]   Telescope, fibers, spectrographs, and CCDs
  - fiberinput[nwave]   Geometric loss at input of fiber for point source.

Required keywords in same HDU EXTNAME=THROUGHPUT (not HDU 0)
  - EXPTIME:  Standard exposure time in seconds
  - EFFAREA:  Effective area of mirror, including obscuration effects, in cm^2
  - FIBERDIA: Fiber diameter in arcsec

Note:
The throughtput format currently does not support per-fiber throughput.

Note:
`fiberinput` in general depends upon source size, seeing, and fiber position
misalignments.  The default is for for perfectly centered point sources.
The format currently does not support extented objects (the user can mimic
this though the input spectra provided to Specter.)

The various forms of throughput affect different types of sources:

                 OBJECT   SKY    CALIB
    EXTINCTION    yes     yes     no
    FIBERINPUT    yes     no      no
    THROUGHPUT    yes     yes     yes

CALIB = calibration lamps or laser systems mounted with the telescope,
i.e. not experiencing sky absorption.


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

Note: Any additional extensions should be read by EXTNAME, not by number;
      the order of other extensions is arbitrary.

HDU 0 keywords:
  - NPIX_X, NPIX_Y : CCD dimensions in pixels

Optional HDU EXTNAME=THROUGHPUT same as THROUGHPUT HDU which could also appear
in a separate fits file.
  
If throughput isn't available, the PSF can still be used to project
photons onto a CCD, but not flux in erg/s/cm^2/A .

Spot Grid PSF
-------------
PSFTYPE = "SPOTGRID"

This PSF type provides an x,y grid of spots to be interpolated.
The first interpolation dimension must monotonically increase with
spectrum number, e.g. the position of a fiber along a slit head.
The second interpolation dimension is wavelength.

HDU 0-2 : Same as Base PSF: X, Y, wavelength or loglam of traces

    HDU SPOTS : spot[i, j, iy, ix]    #- 2D PSF spots
        NAXIS1 = number of spot samples in the spectrum number dimension
        NAXIS2 = number of spot samples in the wavelength direction
        NAXIS3 = size of spot in the CCD y dimension
        NAXIS4 = size of spot in the CCD x direction
    
    HDU SPOTX : spotx[NAXIS1, NAXIS2]   #- CCD X pixel location of spot pixel[0,0]
    HDU SPOTY : spoty[NAXIS1, NAXIS2]   #- CCD Y pixel location of spot pixel[0,0]
    HDU FIBERPOS : fiberpos[nspec]      #- Slit position of each fiber
    HDU SPOTPOS  : spotpos[NAXIS1]      #- Slit positions where spots are sampled
    HDU SPOTWAVE : spotwave[NAXIS2]     #- Wavelengths where spots are sampled

spot[i,j] is a 2D PSF spot sampled at slit position spotpos[i] and
wavelength spotwave[j].  Its center is located on the CCD at
spotx[i,j], spoty[i,j].




























