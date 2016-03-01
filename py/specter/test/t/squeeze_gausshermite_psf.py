#!/usr/bin/env python

"""
Squeeze a DESI GAUSS-HERMITE PSF into something smaller for testing
"""

import sys
import numpy as np
from astropy.io import fits
from astropy.table import Table

infile, outfile = sys.argv

fx = fits.open(infile)

hdr = fx[0].header.copy()
fits.writeto(outfile, np.zeros(0), header=hdr)

hdr = fx[1].header.copy()
nspec = 25
hdr['FIBERMAX'] = nspec-1

data = fx[1].data
tx = Table()
tx['PARAM'] = data['PARAM']
tx['WAVEMIN'] = data['WAVEMIN']
tx['WAVEMAX'] = data['WAVEMAX']
tx['COEFF'] = data['COEFF'][:, 0:nspec, :]

fits.append(outfile, np.array(tx), header=hdr)