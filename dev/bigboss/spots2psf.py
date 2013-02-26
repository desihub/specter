#!/usr/bin/env python

"""
Convert simulated BigBOSS spots into PSF format.

Stephen Bailey, LBL
January 2013
"""

import sys
import os
import numpy as N
from glob import glob
from scipy import ndimage
import fitsio

#+ TODO: refactor this to use numpy.polynomial.legendre instead
from specter.util import LegendreFit

#- Load options
import optparse
parser = optparse.OptionParser(usage = "%prog [options]")
parser.add_option("-i", "--indir", type="string",  help="input directory with spots")
parser.add_option("-o", "--outpsf", type="string",  help="output PSF file")
parser.add_option("-t", "--throughput", type="string",  help="input throughput file to embed with PSF")
parser.add_option("-d", "--debug",  help="start ipython prompt when done", action="store_true")

opts, args = parser.parse_args()

#- for debugging
if opts.indir is None:
    opts.indir = "/data/bigboss/sim/spots/BB_SPEC_20120827difdet/Blue/"

if opts.outpsf is None:
    opts.outpsf = "psf-blat.fits"

#- Hardcode spectrograph and CCD dimensions
CcdPixelSize = 0.015   #- CCD pixel size in mm
SpotPixelSize = 0.001  #- Spot pixel size in mm
FiberSpacing = 0.258   #- center-to-center spacing in mm
GroupSpacing = 0.624   #- center-to-center group gap in mm
FibersPerGroup = 25
GroupsPerSlit = 20
NumFibers = 500
NumPixX = 4096
NumPixY = 4096
nspec = FibersPerGroup * GroupsPerSlit

#- Determine grid of wavelengths and fiber positions for the spots
#- Use set() to get unique values, then convert to sorted array
#- spotgrid maps (fiberpos, wavelength) -> filename
print "Determining wavelength and slit position grid"
wavelength = set()
spotpos = set()
spotfiles = glob(opts.indir + '/*.fits')
spotgrid = dict()
for filename in spotfiles:
    hdr = fitsio.read_header(filename)
    w = hdr['WAVE']*10      #- Wavelength [nm -> AA]
    p = hdr['FIBER']        #- Fiber slit position [mm]
    p = -p      #- Swap slit axis orientation to match CCD x
    wavelength.add(w)       #- Wavelength nm -> AA
    spotpos.add(p)
    spotgrid[(p,w)] = filename
    
wavelength = N.array( sorted(wavelength) )
spotpos = N.array( sorted(spotpos) )

#- Load grid of spots, and the x,y CCD pixel location of those spots
print "Reading spots"
nx = hdr['NAXIS1']
ny = hdr['NAXIS2']
np = len(spotpos)
nw = len(wavelength)
#- spots are mirrored about x=0, thus np -> 2*np + 1
spots = N.zeros( (2*np-1, nw, ny, nx), dtype=N.float32 )
spotx = N.zeros( (np, nw), dtype=N.float32 )
spoty = N.zeros( (np, nw), dtype=N.float32 )
for i, p in enumerate(spotpos):
    for j, w in enumerate(wavelength):
        pix = fitsio.read(spotgrid[(p,w)])
        hdr = fitsio.read_header(spotgrid[(p,w)])        
        
        #- Shift spot to center of image
        #- NOTE: uses spline interpolation, not sinc interpolation
        npy, npx = pix.shape
        yc,xc = ndimage.center_of_mass(pix)
        xmid = (pix.shape[1]-1)/2.0
        ymid = (pix.shape[0]-1)/2.0
        dx = xmid - xc
        dy = ymid - yc
        pix = ndimage.shift(pix, (dy,dx))
                
        spots[i,j] = pix
        if i != len(spotpos)-1:
            spots[-i-1,j] = N.fliplr(pix)
        
        #- Reference pixel in FITS file
        xref = hdr['CRPIX1']-1
        yref = hdr['CRPIX2']-1
        
        #- Location of centroid on CCD in mm from center
        spotx[i,j] = hdr['CRVAL1'] + (xmid-xref+dx)*hdr['CDELT1']
        spoty[i,j] = hdr['CRVAL2'] + (ymid-yref+dy)*hdr['CDELT2']

#- mirror dimensions
blat = N.zeros( (12+11, 23), dtype=spotx.dtype)
blat[0:12] = spotx
blat[12:] = -spotx[-2::-1]
foo = N.zeros( (12+11, 23), dtype=spoty.dtype)
foo[0:12] = spoty
foo[12:] = spoty[-2::-1]

spotx = blat
spoty = foo

#- Convert spotx, spoty to pixel units instead of mm
spotx = spotx/CcdPixelSize + NumPixX/2
spoty = spoty/CcdPixelSize + NumPixY/2

#- Extend spotpos
spotpos = N.concatenate( (spotpos, -spotpos[-2::-1]))
np = len(spotpos)

#- Map location of each fiber along the slit
ifiber = N.arange(NumFibers).astype(int)
ngaps = ifiber / FibersPerGroup    #- Number of gaps prior to fiber ifiber
fiberpos = ifiber*FiberSpacing + ngaps*(GroupSpacing - FiberSpacing)
fiberpos -= N.mean(fiberpos)

#-----
#- Determine range of wavelengths to fit
#- Fit Legendre polynomials and extrapolate to CCD edges
wmin = wavelength[0]
wmax = wavelength[-1]
for i in range(np):
    poly = LegendreFit(spoty[i], wavelength, order=5, xmin=0, xmax=NumPixY)
    wmin = min(wmin, poly(0))
    wmax = max(wmax, poly(NumPixY-1))
    print i, wmin, wmax, poly(0), poly(NumPixY-1)
    
#- Round down/up to nearest 10 AA
wmin = int(wmin/10)*10
wmax = int(wmax/10+1)*10

#-----
#- Determine x,y location of each spectral trace along a wavelength grid
#- `wavelength` is where the spots are sampled; `ww` will be the finer sampling
print "Interpolating X,Y location of traces"
wstep = 2.0
ww = N.arange(wmin, wmax+wstep/2, wstep)

pmin = min(spotpos[0], fiberpos[0])
pmax = max(spotpos[-1], fiberpos[-1])

#- For slices in slit position, fit y vs. w
y_vs_w = N.zeros( (np, len(ww)) )
for i in range(np):
    poly = LegendreFit(wavelength, spoty[i], order=7, xmin=wmin, xmax=wmax)
    y_vs_w[i] = poly(ww)

Y = N.zeros( (nspec, len(ww)), dtype=N.float32 )
for i in range(len(ww)):
    poly = LegendreFit(spotpos, y_vs_w[:, i], order=7, xmin=pmin, xmax=pmax)
    Y[:, i] = poly(fiberpos)

#- for a slice in wavelength, fit x vs. slit position
x_vs_p = N.zeros( (nw, len(fiberpos)) )
for i in range(nw):
    poly = LegendreFit(spotpos, spotx[:,i], order=7, xmin=pmin, xmax=pmax)
    x_vs_p[i] = poly(fiberpos)
    
X = N.zeros( (nspec, len(ww)), dtype=N.float32 )
for i in range(nspec):
    poly = LegendreFit(wavelength, x_vs_p[:, i], order=7, xmin=wmin, xmax=wmax)
    X[i, :] = poly(ww)

#- Wavelength grid
W = N.tile(ww, nspec).reshape( (nspec, len(ww)) ).astype(N.float32)

#-------------------------------------------------------------------------
#- Write to fits file
print "Writing", opts.outpsf

#- Use first spot file for representative header to pass keywords through
hdr = fitsio.read_header(spotfiles[0])
hdr.delete('WAVE')
hdr.delete('FIBER')
hdr.add_record({"name":"PSFTYPE",  "value":"SPOTGRID", "comment":"Grid of simulated PSF spots"})
hdr.add_record({"name":"NPIX_X",   "value":NumPixX,    "comment":"Number of CCD pixels in X direction"})
hdr.add_record({"name":"NPIX_Y",   "value":NumPixY,    "comment":"Number of CCD pixels in Y direction"})
hdr.add_record({"name":"NSPEC",    "value":nspec,      "comment":"Number of spectra"})
hdr.add_record({"name":"NWAVE",    "value":nw,         "comment":"Number of wavelength samples"})
hdr.add_record({"name":"CCDPIXSZ", "value":CcdPixelSize, "comment":"CCD pixel size"})
hdr.add_record({"name":"DFIBER",   "value":FiberSpacing, "comment":"Center-to-center pitch of fibers on slit"})
hdr.add_record({"name":"DGROUP",   "value":GroupSpacing, "comment":"Center-to-center spacing between fiber groups on slit"})
hdr.add_record({"name":"NGROUPS",  "value":GroupsPerSlit,   "comment":"Number of fiber groups per slit"})
hdr.add_record({"name":"NFIBGRP",  "value":FibersPerGroup,  "comment":"Number of fibers per group"})

fitsio.write(opts.outpsf, X, extname='X', header=hdr, clobber=True)

fitsio.write(opts.outpsf, Y, extname='Y')
fitsio.write(opts.outpsf, W, extname='WAVELENGTH')
fitsio.write(opts.outpsf, spots, extname='SPOTS')
fitsio.write(opts.outpsf, spotx, extname='SPOTX')
fitsio.write(opts.outpsf, spoty, extname='SPOTY')
fitsio.write(opts.outpsf, fiberpos, extname='FIBERPOS')
fitsio.write(opts.outpsf, spotpos, extname='SPOTPOS')
fitsio.write(opts.outpsf, wavelength, extname='SPOTWAVE')

#- Add pre-computed throughput to PSF if requested
if opts.throughput:
    header = fitsio.read_header(opts.throughput, 'THROUGHPUT')
    data = fitsio.read(opts.throughput, 'THROUGHPUT')
    fitsio.write(opts.outpsf, data, header=header, extname='THROUGHPUT')

#--- DEBUG ---
if opts.debug:
    import pylab as P
    P.ion()
    import IPython
    IPython.embed()
#--- DEBUG ---





















