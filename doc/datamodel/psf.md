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



