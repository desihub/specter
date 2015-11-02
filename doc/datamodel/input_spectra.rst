====================
Input Spectra Format
====================

Input spectra in FITS format are supported as either image or binary
table HDUs.  The wavelength grid, object type, and flux units are
specified by additional columns, images, and keywords.  Details below.

Input Spectra: Binary Table
---------------------------

HDU 0 : ignored
~~~~~~~~~~~~~~~

HDU 1 : Binary table of spectra
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Required:

* ``flux[nwave]`` or ``flux[nspec, nwave]`` column
* Wavelength grid:
  - header keywords CRVAL1 and CDELT1, or
  - ``wave``, ``wavelength`` or ``loglam`` column with dimensions matching ``flux``

Optional:

* ``objtype[nspec]`` column, or header keyword OBJTYPE.  Default 'STAR'. See below for details on OBJTYPE.
* ``FLUXUNIT`` header keyword (see below)

Other columns and HDUs may be present and will be ignored.

Input Spectra: Image HDU
------------------------

HDU 0 image
~~~~~~~~~~~

* ``flux[nspec, nflux]``
* ``FLUXUNIT`` header keyword (see below)

Wavelength grid defined by one of these:

* HDU 0 keywords ``CRVAL1`` and ``CDELT1`` with optional ``LOGLAM`` flag (see below)
* HDU ``EXTNAME='WAVELENGTH'`` image with wavelength grid in Angstroms
  - ``wavelength[nflux]`` or ``wavelength[nspec, nflux]``
* HDU ``EXTNAME='LOGLAM'`` image with wavelength grid in log10(Angstroms)
  - ``loglam[nflux]`` or ``loglam[nspec, nflux]``

Object type defined by one of these:

* HDU 0 keyword `OBJTYPE` (if all objects are the same)
* HDU `EXTNAME='TARGETINFO'` binary table
  - `OBJTYPE` column required
  - other columns are user-specific and will be ignored

Flux Units
----------

Flux unts are specified by one of the following:

* ``TUNITnn`` keyword for binary table ``flux`` column
* ``FLUXUNIT`` keyword in same HDU as image/table with ``flux``

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
an image HDU, it may be specified by header keywords ``CRVAL1`` and ``CDELT1``
and optionally ``LOGLAM`` (0/1 for linear/log10, default=linear=0).

*e.g.* to specify wavelengths [3600, 3601, 3602, ...]::

    CRVAL1  =                 3600 / Reference wavelength in Angstroms
    CDELT1  =                    1 / delta wavelength

or to specify log10(wavelength) spacing [3.5500, 3.5501, 3.5502, ...]::

    CRVAL1  =               3.5500 / Reference log10(Angstroms)
    CDELT1  =               0.0001 / delta log10(Angstroms)
    LOGLAM  =                    1 / Log10 spaced wavelengths

The FITS standard also requires ``CRPIX1`` and ``CTYPE1`` but these are ignored.

Object Type
-----------

The default object type is "STAR" (astronomical point source), but
other options may be given by either

* ``OBJTYPE`` header keyword, or
* ``objtype[nspec]`` column of binary table format

Options:

* SKY : Geometric losses of input fiber size are not applied
* CALIB : Calibration sources such as arc or flat lamps, or tuneable
  laser systems.  Atmospheric extinction is not applied, nor are
  geometric losses of input fiber size.
* STAR or anything else: treated as an astromical object with all
  sources of throughput loss applied.  
