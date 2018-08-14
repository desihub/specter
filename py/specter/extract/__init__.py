"""
specter.extract
===============

Tools for extracting spectra from 2D images.

Notes
-----

Initial work for adding extractions.  Starting with row-by-row.
This may be renamed and/or refactored into classes::

    from specter.psf import load_psf
    from specter.extract.ex1d import extract1d

    psf = load_psf('data/boss/pixpsf-r1-00140299.fits')
    img = fitsio.read('data/boss/img-r1-00140299.fits', 0)
    ivar = fitsio.read('data/boss/img-r1-00140299.fits', 1)

    specflux, specivar = extract1d(img, ivar, psf, specrange=(20,60), yrange=(1000,1100))

    %time specflux, specivar = extract1d(img, ivar, psf, yrange=(1000,1100), specrange=(20,40) )

==== =====
rows time
==== =====
10   0.32
20   0.59
50   1.75
100  4.94
200  30.7
==== =====

0.59 * 25 * (4000/100) / 60 = 10 minutes.  Not bad. Up to 15 minutes now after refactoring xsigma.

::

    from specter.psf import load_psf
    psf = load_psf('data/boss/pixpsf-r1-00140299.fits')
    psf.xsigma(0, 7000)

    ### yy = N.linspace(10,4000)
    yy = N.arange(10,4000,5)
    ww = psf.wavelength(0, y=yy)
    xsig = [psf.xsigma(0, wavelength=w) for w in ww]

    plot(yy, xsig)

    from specter.util import sincshift
    from scipy.ndimage import center_of_mass
    def xsigma(spot):
        yc, xc = center_of_mass(spot)
        xx = N.arange(spot.shape[1])
        xspot = spot.sum(axis=0)
        return N.sqrt(N.sum(xspot*(xx-xc)**2) / N.sum(xspot))



::

    xr, yr, spot = psf.xypix(0, psf.wavelength(0, y=1000))

    xx = N.arange(xr.start, xr.stop)
    yy = N.arange(yr.start, yr.stop)
    xspot = spot.sum(axis=0)
    yspot = spot.sum(axis=1)

    N.sum(xx*xspot) / N.sum(xspot)
    N.sum(yy*yspot) / N.sum(yspot)

    #- Variance
    from scipy.ndimage import center_of_mass
    yc, xc = center_of_mass(spot)
    xx = N.arange(spot.shape[1])
    yy = N.arange(spot.shape[0])
    N.sum(xspot*(xx-xc)**2) / N.sum(xspot)
    N.sum(yspot*(yy-yc)**2) / N.sum(yspot)
"""

from __future__ import absolute_import, division, print_function, unicode_literals

### from ex1d import ex1d
from .ex2d import ex2d, ex2d_patch
