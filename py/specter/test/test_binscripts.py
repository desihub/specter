#!/usr/bin/env python

"""
Unit tests for executable scripts in specter/bin
"""
from __future__ import print_function
import os
import sys
import numpy as np
import unittest
from astropy.io import fits
from uuid import uuid4

from specter.test import test_data_dir
from specter.psf import load_psf
from specter.extract.ex2d import ex2d
import specter.io

from astropy.io import fits

_base = uuid4().hex
imgfile = 'testimg-'+_base+'.fits'
specfile = 'testspec-'+_base+'.fits'
specfile2 = 'testspec2-'+_base+'.fits'

class TestBinScripts(unittest.TestCase):

    def setUp(self):
        self.specter_dir = os.path.dirname( # top-level
            os.path.dirname( # py/
                os.path.dirname( # specter/
                    os.path.dirname(os.path.abspath(__file__)) # test/
                    )
                )
            )
        self.exspec_cmd = """{executable} {specter_dir}/bin/exspec \
          -i {imgfile} \
          -p {specter_dir}/data/test/psf-monospot.fits \
          -o {specfile} \
          -w 7500,7620,{dwave} \
          --specrange {specmin},{specmax}"""

        #- Add this package to PYTHONPATH so that binscripts can find it
        try:
            self.origPath = os.environ['PYTHONPATH']
            os.environ['PYTHONPATH'] = os.path.join(self.specter_dir,'py') + ':' + self.origPath
        except KeyError:
            self.origPath = None
            os.environ['PYTHONPATH'] = os.path.join(self.specter_dir,'py')
        
    def tearDown(self):
        if self.origPath is None:
            del os.environ['PYTHONPATH']
        else:
            os.environ['PYTHONPATH'] = self.origPath

    def test_aa(self):
        cmd = """{executable} {specter_dir}/bin/specter \
          -i {specter_dir}/data/sky/sky-uves.fits \
          -p {specter_dir}/data/test/psf-monospot.fits \
          -t {specter_dir}/data/test/throughput.fits \
          -o {imgfile} \
          -w 7500,7620 \
          -n -r 0,2 --exptime 1500""".format(
            executable=sys.executable,
            specter_dir=self.specter_dir,
            imgfile = imgfile,
            )
        err = os.system(cmd)
        self.assertEqual(err, 0, 'Error code {} != 0'.format(err))
        self.assertTrue(os.path.exists(imgfile))
        
        fx = fits.open(imgfile)
        self.assertTrue('CCDIMAGE' in fx)
        self.assertTrue('IVAR' in fx)
        fx.close
        
        #- Test the I/O routines while we have the file handy
        image, ivar, hdr = specter.io.read_image(imgfile)
        self.assertTrue(image.shape == ivar.shape)
        
        os.remove(imgfile)
        cmd = cmd + ' --extra'
        err = os.system(cmd)
        self.assertEqual(err, 0, 'Error code {} != 0'.format(err))
        self.assertTrue(os.path.exists(imgfile))
        fx = fits.open(imgfile)
        self.assertTrue('PHOTONS' in fx)
        self.assertTrue('XYWAVE' in fx)
        fx.close

    def test_bb(self):
        for dwave in [1.0, 2.0]:
            cmd = self.exspec_cmd.format(
                executable=sys.executable,
                specter_dir=self.specter_dir,
                imgfile = imgfile,
                specfile = specfile,
                dwave = dwave,
                specmin=0, specmax=2,
                )
            err = os.system(cmd)
            self.assertEqual(err, 0, 'Error code {} != 0 with dwave={}'.format(err, dwave))
            self.assertTrue(os.path.exists(specfile))
            
        fx = fits.open(specfile)
        print(fx.info())
        self.assertTrue('FLUX' in fx)
        self.assertTrue('IVAR' in fx)
        self.assertTrue('WAVELENGTH' in fx)
        self.assertTrue('RESOLUTION' in fx)
        
        #- this is covered in the exspec binscript, but not yet visible to
        #- coverage tools; try it here just for good measure
        specter.io.write_spectra(specfile2,
            fx['WAVELENGTH'].data, fx['FLUX'].data,
            fx['IVAR'].data, fx['RESOLUTION'].data, fx[0].header)
        
        fx.close()
            

    def test_cc(self):
        #- Also check it works for the last fibers and not just the first ones
        cmd = self.exspec_cmd.format(
            executable=sys.executable,
            specter_dir=self.specter_dir,
            imgfile = imgfile,
            specfile = specfile,
            dwave = 1.0,
            specmin=498, specmax=500,
            )
        err = os.system(cmd)
        self.assertEqual(err, 0, 'Error code {} != 0 for --specrange=498,500'.format(err))
        self.assertTrue(os.path.exists(specfile))

    @classmethod
    def tearDownClass(cls):
        for filename in [imgfile, specfile, specfile2]:
            if os.path.exists(filename):
                print("Removing", filename)
                os.remove(filename)


if __name__ == '__main__':
    unittest.main()
