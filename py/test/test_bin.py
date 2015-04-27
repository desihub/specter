#!/usr/bin/env python

"""
Unit tests for executable scripts in specter/bin
"""

import os
import numpy as np
import unittest
from astropy.io import fits
from uuid import uuid4

from specter.test import test_data_dir
from specter.psf import load_psf
from specter.extract.ex2d import ex2d

_base = uuid4().hex
imgfile = 'testimg-'+_base+'.fits'
specfile = 'testspec-'+_base+'.fits'

class TestBin(unittest.TestCase):
        
    def test_aa(self):
        cmd = """specter \
          -i {specter_dir}/data/sky/sky-uves.fits \
          -p {specter_dir}/data/test/psf-monospot.fits \
          -t {specter_dir}/data/test/throughput.fits \
          -o {imgfile} \
          -w 7500,7620 \
          -n -r 0,2 --exptime 1500""".format(
            specter_dir=os.getenv('SPECTER_DIR'),
            imgfile = imgfile,
            )
        err = os.system(cmd)
        self.assertEqual(err, 0, 'Error code {} != 0'.format(err))
        self.assertTrue(os.path.exists(imgfile))

    def test_bb(self):
        cmd = """exspec \
          -i {imgfile} \
          -p {specter_dir}/data/test/psf-monospot.fits \
          -o {specfile} \
          -w 7500,7620,1.0 \
          --specrange 0,2""".format(
            specter_dir=os.getenv('SPECTER_DIR'),
            imgfile = imgfile,
            specfile = specfile,
            )
        err = os.system(cmd)
        self.assertEqual(err, 0, 'Error code {} != 0'.format(err))
        self.assertTrue(os.path.exists(specfile))
        
    @classmethod
    def tearDownClass(cls):
        for filename in [imgfile, specfile]:
            if os.path.exists(filename):
                print "Removing", filename
                os.remove(filename)        
        
            
if __name__ == '__main__':
    unittest.main()
        
