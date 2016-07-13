#!/usr/bin/env python

"""
Test Specter file formats.  Loop over example files and just make sure
that we can read them.
"""

from __future__ import absolute_import, division, print_function, unicode_literals

import os
from os.path import basename, join
from glob import glob
import unittest
from pkg_resources import resource_filename
from ..io import read_simspec

class TestSpecIO(unittest.TestCase):

    def setUp(self):
        self.test_data_dir = resource_filename('specter.test', 't')
        self.specfiles = sorted(glob(join(self.test_data_dir, 'spec-*.fits')))

    def test_files(self):
        wipeout = None
        for specfile in self.specfiles:
            try:
                x = read_simspec(specfile)
            except Exception as e:
                print("Failed on {0}: {1}".format(basename(specfile), str(e)))
                wipeout = e
        if wipeout:
            raise wipeout

if __name__ == '__main__':
    unittest.main()
    # suite = unittest.TestLoader().loadTestsFromTestCase(TestSpecIO)
    # unittest.TextTestRunner(verbosity=2).run(suite)
