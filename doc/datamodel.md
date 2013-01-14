Specter data model
==================

Input spectra
=============

### HDU 0 : blank ###

### HDU 1 : binary table of spectra ###

Required:
  - flux[nwave] or flux[nspec, nwave], or
  - loglam or wavelength (1D[nwave] or 2D[nspec, nwave])
  - keyword FLUXUNIT
    - photons
    - photons/A
    - ergs/s/cm^2/A
    - ergs/s/cm^2/A/arcsec^2

"photons: are treated as photons as observed by the CCD, i.e. all
efficiencies are already applied.  "flux" is treated as "ergs/s/cm"

Optional:
  Header keyword OBJTYPE, or column objtype[nspec]
    values = CALIB, SKY, other (GAL, QSO, STAR)...
    Anything other than CALIB or SKY will be treated as an astronomical source

Other user-specific columns may exist and will be ignored by Specter.

Current implementation: flux is in units of total photons on CCD

Units to support (in progress)
    photons                 - treated as a delta function in wavelength
    photons/A               - to be integrated over a wavelength range
    ergs/s/cm^2/A           - astronomical objects
    ergs/s/cm^2/A/arcsec^2  - sky spectra, calibration lamps/lasers

photons = (flux * EXPTIME * EFFAREA * wavelength / (h*c)) * throughput


Throughput
==========
In Progress

HDU EXTNAME=THROUGHPUT, binary table with columns:
  - wavelength[nwave]   Wavelength in Angstroms
      or loglam[nwave]    or log10(Angstroms)
  - extinction[nwave]   Atmospheric extinction in mags/airmass
                        Atm throughput = 10^(0.4 * extinction * airmass)
  - throughput[nwave]   Telescope, fibers, spectrographs, and CCDs
  - fiberinput[nwave]   Geometric loss at input of fiber for point source.

TODO: throughput and fiberinput may optinally have dimensions [nfiber, nwave]

Keywords in same HDU EXTNAME=THROUGHPUT
  - EXPTIME:  Standard exposure time in seconds
  - EFFAREA:  Effective area of mirror, including obscuration effects, in cm^2
  - FIBERDIA: Fiber diameter in arcsec

  flux units = ergs/s/cm^2/A
  photons/A = flux * EXPTIME * EFFAREA * wavelength / (h*c)

fiberinput in general depends upon source size, seeing, and fiber position
misalignments.  Default is for for perfectly centered point source.
--> This is the weak point of this format.

The various forms of throughput affect different types of sources:

                 Object   Sky    Calib
    EXTINCTION    yes     yes     no
    FIBERINPUT    yes     no      no
    THROUGHPUT    yes     yes     yes

Calib = calibration lamps or laser systems.



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

--> Additional extensions should be read by EXTNAME, not by number;
    the order of these extensions is arbitrary

HDU 0 keywords:
  - NPIX_X, NPIX_Y : CCD dimensions in pixels

HDU EXTNAME=THROUGHPUT same as THROUGHPUT HDU which could also appear
in a separate fits file.
  
If throughput isn't available, the PSF can still be used to project
photons onto a CCD, but not flux in ergs/s/cm^2/A .

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




























