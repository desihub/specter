#!/usr/bin/env python

"""
Test specter file formats
"""

import os
from os.path import basename
from glob import glob
import specter.io
import unittest

class TestSpecIO(unittest.TestCase):
    def setUp(self):
        indir = os.environ['SPECTER_DIR'] + '/test/data'
        self.specfiles = sorted(glob(indir+'/spec-*.fits'))
        
    def test_files(self):
        wipeout = None
        for specfile in self.specfiles:
            try:
                x = specter.io.read_simspec(specfile)
            except Exception, e:
                print "Failed on %s: %s" % (basename(specfile), str(e))
                wipeout = e
        if wipeout:
            raise wipeout
            
if __name__ == '__main__':
    unittest.main()            
    # suite = unittest.TestLoader().loadTestsFromTestCase(TestSpecIO)
    # unittest.TextTestRunner(verbosity=2).run(suite)
        