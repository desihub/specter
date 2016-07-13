#!/usr/bin/env python

"""
Take a BigBOSS spotgrid PSF and trim it down for testing
"""

from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import os
import numpy as np
from astropy.io import fits

infile = os.environ['SPECTER_DIR'] + "/data/bigboss/designs/20120827difdet/bbpsf-I.fits"
outfile = os.environ['SPECTER_DIR'] + "/test/data/psf-spot.fits"

fx = fits.open(infile)
x, xhdr = fx[0].data, fx[0].header
y, yhdr = fx[1].data, fx[1].header
w, whdr = fx[2].data, fx[2].header

spots, spotshdr = fx[3].data, fx[3].header
spotx, spotxhdr = fx[4].data, fx[4].header
spoty, spotyhdr = fx[5].data, fx[5].header

fiberpos, fposhdr   = fx[6].data, fx[6].header
spotpos, sposhdr    = fx[7].data, fx[7].header
spotwave, swavehdr  = fx[8].data, fx[8].header
throughput, thruhdr = fx[9].data, fx[9].header

fits.writeto(outfile, x[:, 0::10], header=xhdr, clobber=True)
fits.append(outfile, y[:, 0::10], header=yhdr)
fits.append(outfile, w[:, 0::10], header=whdr)

fits.append(outfile, spots[0::5, 0::5, :, :], header=spotshdr)
fits.append(outfile, spotx[0::5, 0::5], header=spotxhdr)
fits.append(outfile, spoty[0::5, 0::5], header=spotyhdr)

fits.append(outfile, fiberpos, header=fposhdr)
fits.append(outfile, spotpos[0::5], header=sposhdr)
fits.append(outfile, spotwave[0::5], header=swavehdr)

fits.append(outfile, throughput, header=thruhdr)

