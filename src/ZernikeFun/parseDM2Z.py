import numpy
import scipy
import pyfits

data = pyfits.open('/home/deen/Data/GRAVITY/SPARTA_Data/Casey/Matrices/DM2Z.fits')

newData = numpy.array(data[0].data[0:50,:], dtype=numpy.float32)

HDU = pyfits.PrimaryHDU(newData)
HDU.writeto('LoopMonitor.DMPOS2MODES.fits', clobber=True)
