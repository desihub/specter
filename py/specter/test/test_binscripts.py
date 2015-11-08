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

_base = uuid4().hex
imgfile = 'testimg-'+_base+'.fits'
specfile = 'testspec-'+_base+'.fits'

class TestBinScripts(unittest.TestCase):

    def setUp(self):
        self.specter_dir = os.path.dirname( # top-level
            os.path.dirname( # py/
                os.path.dirname( # specter/
                    os.path.dirname(__file__) # test/
                    )
                )
            )
        self.exspec_cmd = """{executable} {specter_dir}/bin/exspec \
          -i {imgfile} \
          -p {specter_dir}/data/test/psf-monospot.fits \
          -o {specfile} \
          -w 7500,7620,{dwave} \
          --specrange {specmin},{specmax}"""

    @unittest.skip("Scripts need to be refactored for test purposes.")
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

    @unittest.skip("Scripts need to be refactored for test purposes.")
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

    @unittest.skip("Scripts need to be refactored for test purposes.")
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
        for filename in [imgfile, specfile]:
            if os.path.exists(filename):
                print("Removing", filename)
                os.remove(filename)


if __name__ == '__main__':
    unittest.main()
