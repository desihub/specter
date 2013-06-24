Specter Data Model
==================

Input Spectra Format
====================

Input spectra in FITS format are supported as either image or binary
table HDUs.  The wavelength grid, object type, and flux units are
specified by additional columns, images, and keywords.  Details below.

Input Spectra: Binary Table
---------------------------

### HDU 0 : ignored ###

### HDU 1 : Binary table of spectra ###

Required:

  * `flux[nwave]` or `flux[nspec, nwave]` column
  * Wavelength grid:
    - header keywords CRVAL1 and CDELT1, or
    - `wave`, `wavelength` or `loglam` column with dimensions matching `flux`

Optional:

  * `objtype[nspec]` column, or header keyword OBJTYPE.  Default 'STAR'
    See below for details on OBJTYPE.
  * `FLUXUNIT` header keyword (see below)
  
Other columns and HDUs may be present and will be ignored.


Input Spectra: Image HDU
------------------------

HDU 0 image : `flux[nspec, nflux]`

  * `FLUXUNIT` header keyword (see below)

Wavelength grid defined by one of these:

  * HDU 0 keywords `CRVAL1` and `CDELT1` with optional `LOGLAM` flag
    (see below)
  * HDU `EXTNAME='WAVELENGTH'` image with wavelength grid in Angstroms
    - `wavelength[nflux]` or `wavelength[nspec, nflux]`
  * HDU `EXTNAME='LOGLAM'` image with wavelength grid in log10(Angstroms)
    - `loglam[nflux]` or `loglam[nspec, nflux]`
  
Object type defined by one of these:

  * HDU 0 keyword `OBJTYPE` (if all objects are the same)
  * HDU `EXTNAME='TARGETINFO'` binary table
    - `OBJTYPE` column required
    - other columns are user-specific and will be ignored

Flux Units
----------

Flux unts are specified by one of the following:

  * `TUNITnn` keyword for binary table `flux` column
  * `FLUXUNIT` keyword in same HDU as image/table with `flux`

Options are:

  * Treated as function values to be multipled by bin width:
    - erg/s/cm^2/A  (default)
    - erg/s/cm^2/A/arcsec^2
    - photon/A
    - photon/A/arcsec^2
  * Treated as delta functions at each given wavelength:
    - photons
    - erg/s/cm^2
    - erg/s/cm^2/arcsec^2

For example, an astromical object is typically in units "erg/s/cm^2/A"
and will be converted to photons using all throughput terms of the
throughput model.  A sky spectrum may be in "erg/s/cm^2/A/arcsec^2" and
will be multiplied by the area of the fiber instead of having a
fiber input geometric loss applied.

Wavelength Grid
---------------

If the wavelength grid is not specified by a binary table column or
an image HDU, it may be specified by header keywords `CRVAL1` and `CDELT1`
and optionally `LOGLAM` (0/1 for linear/log10, default=linear=0).

e.g. to specify wavelengths [3600, 3601, 3602, ...]:

    CRVAL1  =                 3600 / Reference wavelength in Angstroms
    CDELT1  =                    1 / delta wavelength

or to specify log10(wavelength) spacing [3.5500, 3.5501, 3.5502, ...]

    CRVAL1  =               3.5500 / Reference log10(Angstroms)
    CDELT1  =               0.0001 / delta log10(Angstroms)
    LOGLAM  =                    1 / Log10 spaced wavelengths

The FITS standard also requires `CRPIX1` and `CTYPE1` but these are ignored.

Object Type
-----------

The default object type is "STAR" (astronomical point source), but
other options may be given by either

  * `OBJTYPE` header keyword, or
  * `objtype[nspec]` column of binary table format

Options:

  * SKY : Geometric losses of input fiber size are not applied
  * CALIB : Calibration sources such as arc or flat lamps, or tuneable
    laser systems.  Atmospheric extinction is not applied, nor are
    geometric losses of input fiber size.
  * STAR or anything else: treated as an astromical object with all
    sources of throughput loss applied.  


Output Image Format
===================

### HDU 0 CCDIMAGE ###

Image of spectra projected onto the CCD with the PSF, with optional noise.
The readout noise is always Gaussian, but the photon noise could be
Poisson (default) or Gaussian (if --gaussnoise option was used.)

Header keywords:

  - SIMDATA = True
  - PREPROC = True
  - GAIN    = CCD gain in electrons/ADU
  - RDNOISE = CCD amplifier readout noise in electrons
  - SIMNOISE = "Gaussian", "Poisson", or "None"

### HDU 1 IVAR ###

If noise was added to the CCD, this contains the pixel inverse variance
if the noise was treated as Gaussian (even if the photon shot noise is
Poisson).

### HDU 2 PHOTONS ###

*NOTE* : The format of this HDU will likely change to match the format
of input spectra (`flux` plus `FLUXUNIT` keyword instead of `PHOTONS`).
This has not yet been implemented.

Optional HDU with binary table giving the photon spectra projected onto
the CCD after all throughputs were applied.  Extraction code should
produce this before any calibrations.  There is one row per spectrum,
with columns:

  - PHOTONS[nwave]
  - WAVELENGTH[nwave]

### HDU 3 XYWAVE ###

Optional HDU with x,y vs. wavelength of spectral traces on CCD.
There is one row per spectrum, with columns:

  - X[nwave]
  - Y[nwave]
  - WAVELENGTH[nwave]


Throughput Format
=================

### HDU EXTNAME=THROUGHPUT ###

binary table with columns:

  - `wavelength[nwave]`   Wavelength in Angstroms,
      or `loglam[nwave]`  for log10(Angstroms)
  - `extinction[nwave]`   Atmospheric extinction in mags/airmass;
                          Atm throughput = 10^(0.4 * extinction * airmass)
  - `throughput[nwave]`   Telescope, fibers, spectrographs, and CCDs
  - `fiberinput[nwave]`   Geometric loss at input of fiber for point source.

Required keywords in same HDU `EXTNAME=THROUGHPUT` (not HDU 0)

  - `EXPTIME`:  Standard exposure time in seconds
  - `EFFAREA`:  Effective area of mirror, including obscuration effects, in cm^2
  - `FIBERDIA`: Fiber diameter in arcsec

*Note*:
The throughtput format currently does not support per-fiber throughput.

*Note*:
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

  - NPIX\_X, NPIX\_Y : CCD dimensions in pixels

Optional HDU EXTNAME=THROUGHPUT in the same format as described above.
i.e. the throughput may be kept in a separate FITS file, or bundled with the
instrument PSF for convenience.
  
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

`spot[i,j]` is a 2D PSF spot sampled at slit position `spotpos[i]` and
wavelength `spotwave[j]`.  Its center is located on the CCD at
`spotx[i,j], spoty[i,j]`.

Pixellated PSF PCA expansion
----------------------------
PSFTYPE = "PCA-PIX"

This format is a PCA-like pixelated model such that

    pix = ConstImage + x*ImageX + y*ImageY + x*y*ImageXY + ...

HDU 0-2 : Same as Base PSF: X, Y, wavelength or loglam of traces

HDU 3 : table with x and y exponents to use for each image model

    Columns IMODEL, XEXP, YEXP
    One row per model image

HDU 4 : table with x and y scale factors

    Columns IGROUP, X0, XSCALE, Y0, YSCALE
    One row per fiber
    
The model images are separately calculated for groups of fibers
(could be bundles, but doesn't have to be).  Within each group,
the x and y dimensions are rescaled such that

    x = xscale * (xpix - x0)
    y = yscale * (ypix - y0)

HDU 5 : 4-D array with the PSF images

    Array[IGROUP, IMODEL, IY, IX]
    
e.g. to find the PSF for ifiber iflux:

    xpix = HDU0[ifiber, iflux]
    ypix = HDU1[ifiber, iflux]
    x0 = HDU4.X0[ifiber]
    y0 = HDU4.Y0[ifiber]
    xscale = HDU4.XSCALE[ifiber]
    yscale = HDU4.YSCALE[ifiber]    
    x = xscale*(xpix - x0)
    y = yscale*(ypix - y0)
    
    igroup = HDU4.IGROUP[ifiber]
    psf = [BlankImage]
    for imodel in range( len(HDU3) ):
        xexp = HDU3.XEXP[imodel]
        yexp = HDU3.YEXP[imodel]
        psf += x^xexp * y^yexp * HDU5[igroup, imodel]

Other PSF Formats
-----------------

Specter grew out of "bbspec" which includes PSF formats for:

  * 2D rotated asymmetric Gaussian
  * Gauss-Hermite polynomials

These are structurally compatible with Specter PSFs but haven't been
ported yet.




























