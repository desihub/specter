#!/usr/bin/env python

"""
Convert templates from BOSS spZbest files into Specter input spectra files.
The resolution won't be correct, but this can provide a useful set of
templates of a diverse class of targets all mixed together in a single file.

Stephen Bailey
January 2013
"""

import sys
import os
import numpy as N
import fitsio

import optparse

parser = optparse.OptionParser(usage = "%prog [options]")
parser.add_option("-i", "--input", type="string",  help="input spZbest")
parser.add_option("-o", "--output", type="string",  help="output file")
# parser.add_option("-x", "--xxx",   help="some flag", action="store_true")

opts, args = parser.parse_args()

#- Read data
hdr = fitsio.read_header(opts.input, 0)
zbest = fitsio.read(opts.input, 1).view(N.recarray)
templates = fitsio.read(opts.input, 2)
nspec, nwave = templates.shape

#- Extract columns of interest
z = zbest.Z
zerr = zbest.Z_ERR
zwarn = zbest.ZWARNING
klass = zbest.CLASS

#- Replace OBJTYPE='GALAXY' entries with NOQSO versions of results
igal = (zbest.OBJTYPE == 'GALAXY')
z[igal] = zbest.Z_NOQSO[igal]
zerr[igal] = zbest.Z_ERR_NOQSO[igal]
zwarn[igal] = zbest.ZWARNING_NOQSO[igal]
klass[igal] = zbest.CLASS[igal]

#- Replace SKY and SPECTROPHOTO_STD targets
isky = (zbest.OBJTYPE == "SKY")
klass[isky] = 'SKY'
templates[isky, :] = 0.0

istd = (zbest.OBJTYPE == "SPECTROPHOTO_STD")
klass[istd] = 'STAR'

#- Create output array
a = N.array(zip(templates, z, zerr, zwarn, klass),
            dtype=[('FLUX',     str(templates.dtype), (nwave,)),
                   ('Z',        str(z.dtype)),
                   ('ZERR',     str(zerr.dtype)),
                   ('ZWARNING', str(zwarn.dtype)),
                   ('OBJTYPE',  str(klass.dtype)),
                   ])

#- Write output table                   
fitsio.write(opts.output, a, clobber=True)

#- Update keywords
fx = fitsio.FITS(opts.output, 1)
fx[1].write_key('CRVAL1', hdr['CRVAL1'], 'Starting log10(Angstrom) value')
fx[1].write_key('CDELT1', hdr['CD1_1'],  'delta log10(Angstrom) step')
fx[1].write_key('DC-FLAG', 1,  'log10(Angstrom) scale')
fx[1].write_key('TUNIT1', "1e-17 erg/s/cm^2/A",  'Flux units')
fx.close()





