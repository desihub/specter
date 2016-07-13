#!/usr/bin/env python

"""
1D Extraction like Horne 1986

Stephen Bailey, LBL
Spring 2013
"""

from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import os
import numpy as np
import math

from specter.util import gausspix, weighted_solve

def ex1d(img, mask, psf, readnoise=2.5,
              specrange=None, yrange=None,
              nspec_per_group=20, debug=False, model=False):
    """
    Extract spectra from an image using row-by-row weighted extraction.
    
    Inputs:
        img[ny, nx]     CCD image
        mask[ny, nx]    0=good, non-zero=bad
        psf object
        
    Optional Inputs:
        readnoise = CCD readnoise
        specrange = (specmin, specmax) Spectral range to extract (default all)
        yrange = (ymin, ymax) CCD y (row) range to extract (default all rows)
        --> ranges are python-like, i.e. yrange=(0,100) extracts 100 rows
            from 0 to 99 inclusive but not row 100.
            
        groupspec: extract spectra in groups of N spectra
            (faster if spectra are physically separated into non-overlapping
            groups)
            
        debug: if True, stop with prompt after each row
        
    Returns:
        spectra[nspec, ny]   - extracted spectra
        specivar[nspec, ny]  - inverse variance of spectra  
    """

    #- Range of spectra to extract
    specmin, specmax = specrange if (specrange is not None) else (0, psf.nspec)        
    nspec = specmax - specmin
    
    #- Rows to extract
    ymin, ymax = yrange if (yrange is not None) else (0, psf.npix_y)
    ny = ymax - ymin
    nx = img.shape[0]
    xx = np.arange(nx)
    
    spectra = np.zeros((nspec, ny))
    specivar = np.zeros((nspec, ny))
    
    if model:
        imgmodel = np.zeros(img.shape)
        
    #- Loop over groups of spectra
    for speclo in range(specmin, specmax, nspec_per_group):
        spechi = min(specmax, speclo+nspec_per_group)
                
        #- Calc trace centers (x0) and gaussian sigmas (xsigma) for each row
        allx0 = np.zeros((nspec_per_group, ymax-ymin))
        allxsigma = np.zeros((nspec_per_group, ymax-ymin))
        rows = np.arange(ymin, ymax)
        for ispec in range(speclo, spechi):
            w = psf.wavelength(ispec, y=rows)
            allx0[ispec-speclo] = psf.x(ispec, w)
            allxsigma[ispec-speclo] = psf.xsigma(ispec, w)
                
        #- Loop over CCD rows
        for irow, row in enumerate(range(ymin, ymax)):
            if debug and row%500 == 0:
                print("Row {:3d} spectra {}:{}".format(row, speclo, spechi))
        
            #- Determine x range covered for this row of this group of spectra
            wlo = psf.wavelength(speclo, y=row)
            whi = psf.wavelength(spechi-1, y=row)
            if speclo == 0:
                xmin = 0
            else:
                xmin = int(0.5*(psf.x(speclo-1, wlo) + psf.x(speclo, wlo)))
        
            if spechi >= psf.nspec:
                xmax = psf.npix_x
            else:
                xmax = int(0.5*(psf.x(spechi-1, wlo) + psf.x(spechi, wlo)) + 1)
                
            #- Design matrix for pixels = A * flux for this row
            A = np.zeros( (xmax-xmin, spechi-speclo) )
            for ispec in range(speclo, spechi):
                # w = psf.wavelength(ispec, y=row)
                # x0 = psf.x(ispec, w)
                # xsigma = psf.xsigma(ispec, w)
                x0 = allx0[ispec-speclo, irow]
                xsigma = allxsigma[ispec-speclo, irow]
            
                #- x range for single spectrum on single row
                xlo = max(xmin, int(x0-5*xsigma))
                xhi = min(xmax, int(x0+5*xsigma+1))
            
                A[xlo-xmin:xhi-xmin, ispec-speclo] = gausspix(xx[xlo:xhi], x0, xsigma)
            
            #- Original solve, weighting by input ivar
            ### tmpspec, iCov = weighted_solve(A, img[row, xmin:xmax], ivar[row, xmin:xmax])

            #- Solve weighting only by readnoise and mask
            nx = xmax-xmin
            xvar = np.ones(nx)*readnoise**2 * (mask[row, xmin:xmax] == 0)
            tmpspec, iCov = weighted_solve(A, img[row, xmin:xmax], 1.0/xvar)
            
            #- Re-extract with weight incluing model shot noise
            xvar = (A.dot(tmpspec) + readnoise**2) * (mask[row, xmin:xmax] == 0)
            tmpspec, iCov = weighted_solve(A, img[row, xmin:xmax], 1.0/xvar)

            if model:
                imgmodel[row, xmin:xmax] = A.dot(tmpspec)
                
            spectra[speclo-specmin:spechi-specmin, row-ymin] = tmpspec
            specivar[speclo-specmin:spechi-specmin, row-ymin] = iCov.diagonal()   
            
            if debug:
                print(speclo, row)
                import IPython
                IPython.embed()                
                
    if model:
        return spectra, specivar, imgmodel
    else:
        return spectra, specivar
    
        
