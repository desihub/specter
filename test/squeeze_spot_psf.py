#!/usr/bin/env python

"""
Take a BigBOSS spotgrid PSF and trim it down for testing
"""

import sys
import os
import numpy as N
import fitsio

infile = os.environ['SPECTER_DIR'] + "/data/bigboss/designs/20120827difdet/bbpsf-I.fits"
outfile = os.environ['SPECTER_DIR'] + "/test/data/psf-spot.fits"

fx = fitsio.FITS(infile)
x, xhdr = fx[0].read(), fx[0].read_header()
y, yhdr = fx[1].read(), fx[1].read_header()
w, whdr = fx[2].read(), fx[2].read_header()

spots, spotshdr = fx[3].read(), fx[3].read_header()
spotx, spotxhdr = fx[4].read(), fx[4].read_header()
spoty, spotyhdr = fx[5].read(), fx[5].read_header()

fiberpos, fposhdr  = fx[6].read(), fx[6].read_header()
spotpos, sposhdr   = fx[7].read(), fx[7].read_header()
spotwave, swavehdr = fx[8].read(), fx[8].read_header()
throughput, thruhdr = fx[9].read(), fx[9].read_header()

fitsio.write(outfile, x[:, 0::10], header=xhdr, clobber=True)
fitsio.write(outfile, y[:, 0::10], header=yhdr)
fitsio.write(outfile, w[:, 0::10], header=whdr)

fitsio.write(outfile, spots[0::5, 0::5, :, :], header=spotshdr)
fitsio.write(outfile, spotx[0::5, 0::5], header=spotxhdr)
fitsio.write(outfile, spoty[0::5, 0::5], header=spotyhdr)

fitsio.write(outfile, fiberpos, header=fposhdr)
fitsio.write(outfile, spotpos[0::5], header=sposhdr)

