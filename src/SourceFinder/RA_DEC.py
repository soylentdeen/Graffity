import scipy
import numpy
from astropy.coordinates import SkyCoord, ICRS, FK5
import astropy.units as u
import sys

c1 = SkyCoord(sys.argv[1], sys.argv[2], unit='deg')
c2 = SkyCoord(sys.argv[3], sys.argv[4], unit='deg')

print c1.fk5.to_string('hmsdms')
print c2.fk5.to_string('hmsdms')
print c1.separation(c2).arcsec


