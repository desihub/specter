#!/usr/bin/env python

"""
Test specter throughput file format
"""

import os
import numpy as N
import unittest

from specter.throughput import load_throughput
from specter.test import test_data_dir

class TestThroughput(unittest.TestCase):

    def setUp(self):
        self.thru = load_throughput(test_data_dir()+'/throughput.fits')
        self.w = N.arange(5000, 9000, 1)
        self.flux = N.random.uniform(1,2, size=self.w.shape) * 1e-17
        self.photflux = N.random.uniform(0,1, size=self.w.shape)
        
    def test_effarea(self):
        self.assertTrue(self.thru.effarea > 0)

    def test_exptime(self):
        self.assertTrue(self.thru.exptime > 0)

    def test_fiberdia(self):
        self.assertTrue(self.thru.fiberdia > 0)
        
    def test_atmthru(self):
        t = self.thru.atmospheric_throughput(self.w)
        self.assertTrue(N.all( (0.0<=t) & (t<= 1.0)))

    def test_atmext(self):
        ext = self.thru.extinction(self.w)
        self.assertTrue(N.all( ext > 0 ))

    def test_fiberin(self):
        t = self.thru.fiberinput_throughput(self.w)
        self.assertTrue(N.all( (0.0<=t) & (t<= 1.0)))

    def test_hardware(self):
        t = self.thru.hardware_throughput(self.w)
        self.assertTrue(N.all( (0.0<=t) & (t<= 1.0)))
        
    def test_calibthru(self):
        t1 = self.thru(self.w, objtype='CALIB', airmass=1.0)
        t2 = self.thru(self.w, objtype='CALIB', airmass=2.0)
        self.assertTrue( N.all(t1 == t2) )
        self.assertTrue( len(t1) == len(t2) )
        self.assertTrue( len(t1) == len(self.w) )
        self.assertTrue( N.any(t2>0.0) )

    def test_skythru(self):
        t1 = self.thru(self.w, objtype='SKY', airmass=1.0)
        t2 = self.thru(self.w, objtype='SKY', airmass=2.0)
        self.assertTrue( N.all( (t1 > t2) | (t2 == 0.0)) )
        self.assertTrue( len(t1) == len(t2) )
        self.assertTrue( len(t1) == len(self.w) )
        self.assertTrue( N.any(t2>0.0) )

    def test_starthru(self):
        t1 = self.thru(self.w, objtype='STAR', airmass=1.0)
        t2 = self.thru(self.w, objtype='STAR', airmass=2.0)
        self.assertTrue( N.all( (t1 > t2) | (t2 == 0.0)) )
        self.assertTrue( len(t1) == len(t2) )
        self.assertTrue( len(t1) == len(self.w) )
        self.assertTrue( N.any(t2>0.0) )
        
    def test_thru(self):
        t1 = self.thru(self.w, objtype='CALIB', airmass=1.1)
        t2 = self.thru(self.w, objtype='SKY', airmass=1.1)
        t3 = self.thru(self.w, objtype='STAR', airmass=1.1)
        
        self.assertTrue( N.all( (t1>t2) | (t1==0.0) ) )
        self.assertTrue( N.all( (t2>t3) | (t2==0.0) ) )
        self.assertTrue( N.any(t3>0.0) )

    def test_fluxunits(self):
        units = [
            "erg/s/cm^2",
            "erg/s/cm^2/A",
            "erg/s/cm^2/A/arcsec^2",
            "erg/s/cm^2/arcsec^2",
        ]
        for u in units:
            p = self.thru.photons(self.w, self.flux, units=u, objtype='STAR')
            self.assertTrue( N.any(p>0) and N.all(p>=0) )

    def test_scaledunits(self):
        scale = 1e-16
        units = [
            "erg/s/cm^2",
            "erg/s/cm^2/A",
            "erg/s/cm^2/A/arcsec^2",
            "erg/s/cm^2/arcsec^2",
        ]
        for u in units:
            scaled_units = str(scale) + " " + u
            p0 = self.thru.photons(self.w, self.flux, units=u, objtype='STAR')
            p1 = self.thru.photons(self.w, self.flux/scale, units=scaled_units, objtype='STAR')
            ii = (p0 != 0.0)
            dp = N.abs( (p0-p1)[ii]/p0[ii] )
            self.assertTrue( N.max(dp) < 1e-14 )  #- allow for roundoff error

    def test_photunits(self):
        units = [
            "photon",
            "photon/A",
            ### "photon/A/arcsec^2",
        ]
        for u in units:
            p = self.thru.photons(self.w, self.photflux, units=u, objtype='STAR')
            self.assertTrue( N.any(p>0) and N.all(p>=0) )
        
    def test_calibphot(self):
        p1 = self.thru.photons(self.w, self.flux, objtype='CALIB', airmass=1.0)
        p2 = self.thru.photons(self.w, self.flux, objtype='CALIB', airmass=2.0)
        self.assertTrue( N.any(p1>0) and N.all(p1==p2) )

    def test_skyphot(self):
        p1 = self.thru.photons(self.w, self.flux, objtype='SKY', airmass=1.0)
        p2 = self.thru.photons(self.w, self.flux, objtype='SKY', airmass=2.0)
        self.assertTrue( N.any(p1>0) )
        self.assertTrue( N.all( (p1>p2) | (p2==0) ) )

    def test_objphot(self):
        p1 = self.thru.photons(self.w, self.flux, objtype='STAR', airmass=1.0)
        p2 = self.thru.photons(self.w, self.flux, objtype='STAR', airmass=2.0)
        self.assertTrue( N.any(p1>0) )
        self.assertTrue( N.all( (p1>p2) | (p2==0) ) )
        
    def test_multiphot(self):
        p1 = self.thru.photons(self.w, self.flux, objtype='CALIB', airmass=1.0)
        p2 = self.thru.photons(self.w, self.flux, objtype='SKY', airmass=1.0)
        p3 = self.thru.photons(self.w, self.flux, objtype='STAR', airmass=1.0)
        
        self.assertTrue( N.all( (p1>p2) | (p1==0.0) ) )
        self.assertTrue( N.all( (p2>p3) | (p2==0.0) ) )
        self.assertTrue( N.any(p3>0.0) )

    def test_apply_throughput(self):
        f1 = self.thru.apply_throughput(self.w, self.flux, objtype='CALIB')
        f2 = self.thru.apply_throughput(self.w, self.flux, objtype='SKY')
        f3 = self.thru.apply_throughput(self.w, self.flux, objtype='STAR')
    
        self.assertTrue( N.all( (f1>f2) | (f1==0.0) ) )
        self.assertTrue( N.all( (f2>f3) | (f2==0.0) ) )
        self.assertTrue( N.any(f3>0.0) )
        
        self.assertTrue( N.all(f1 <= self.flux) )
        self.assertTrue( N.all(f2 <= self.flux) )
        self.assertTrue( N.all(f3 <= self.flux) )
    
        self.assertTrue( N.all(f2 <= f1) )
        self.assertTrue( N.all(f3 <= f2) )

    #- Throughputs don't work when the model goes negative.
    #- Should we have less negative flux, or even more negative flux?
    #- The can happen even when the integrated flux in a bin is positive
    #- but the underlying model fit goes negative.
    #- This test is a placeholder to remind us of this issue.
    @unittest.expectedFailure
    def test_low_flux(self):
        flux = N.random.uniform(0, 1, len(self.w))
        f1 = self.thru.apply_throughput(self.w, flux, objtype='CALIB')
        f2 = self.thru.apply_throughput(self.w, flux, objtype='SKY')
        f3 = self.thru.apply_throughput(self.w, flux, objtype='STAR')

        self.assertTrue( N.all(f1 <= flux) )
        self.assertTrue( N.all(f2 <= flux) )
        self.assertTrue( N.all(f3 <= flux) )

        self.assertTrue( N.all(f2 <= f1) )
        self.assertTrue( N.all(f3 <= f2) )
            
if __name__ == '__main__':
    unittest.main()            
    # suite = unittest.TestLoader().loadTestsFromTestCase(TestThroughputIO)
    # unittest.TextTestRunner(verbosity=2).run(suite)
        