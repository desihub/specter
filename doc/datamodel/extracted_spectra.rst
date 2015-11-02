========================
Extracted Spectra Format
========================

There are 4 HDUs with EXTNAME='FLUX', 'IVAR', 'WAVELENGTH', and 'RESOLUTION'

HDU EXTNAME='FLUX'
------------------

2D image with ``flux[ispec, wavelength]`` in photons per Angstrom.
No fiber flat fielding or any calibrations have been applied.

Header keywords are propaged from the input image, with the following changes/additions::

    NAXIS1  = Number of wavelength bins
    NAXIS2  = Number of spectra
    SPECMIN = First spectrum
    SPECMAX = Last spectrum (inclusive)
    NSPEC   = NAXIS2 = Number of spectra
    WAVEMIN = First wavelength [Angstroms]
    WAVEMAX = Last wavelength [Angstroms]
    SPECTER = Specter version
    IN_PSF  = Input spectral PSF
    IN_IMG  = Input image

HDU EXTNAME='IVAR'
------------------

Inverse variance of extracted spectra, same dimensions as FLUX.
Units (photons/A)^-2.

HDU EXTNAME='WAVELENGTH'
------------------------

1D array with the wavelength grid in Angstroms.  NAXIS1 same as
'FLUX' and 'IVAR' HDUs.

HDU EXTNAME='RESOLUTION'
------------------------

3D array storing the per-fiber band-diagonal elements of the spectroperfectionism "Resolution Matrix"::

    NAXIS1 = Number of wavelength bins
    NAXIS2 = Number of diagonal bands
    NAXIS3 = Number of spectra

The format is designed such that the following python code will create a sparse matrix representation of the resolution matrix R for spectrum i::

    nwave = header['NAXIS1']
    nband = header['NAXIS2']
    offsets = range(nband//2,-nband//2,-1)
    R = scipy.sparse.dia_matrix((data[i], offsets), (nwave,nwave))

Also see :download:`ex2d_ResMatrix.pdf` for how this is created from the divide-and-conquer extractions.

Cross terms between fibers are not stored.
