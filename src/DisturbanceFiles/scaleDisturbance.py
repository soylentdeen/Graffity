import scipy
import numpy
import pyfits
import sys

data = pyfits.getdata('/home/deen/Data/GRAVITY/SPARTA_Data/DisturbanceFiles/cimdatDMNoise.fits')

newData = data/4.0

newHDU = pyfits.PrimaryHDU(data=numpy.array(newData, dtype=numpy.float32))

newHDU.writeto('/home/deen/Data/GRAVITY/SPARTA_Data/DisturbanceFiles/newDisturb.fits', clobber=True)

