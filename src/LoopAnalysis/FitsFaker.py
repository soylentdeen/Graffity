import pyfits
import numpy
import scipy

angleStart = 0.0
angleStop = 180.0
deltaAngle = 1.0

data = numpy.arange(angleStart, angleStop, deltaAngle, dtype=numpy.float32)

outfile = 'RTC.ROTLIB.ANGLE_TABLE.fits'

hdu = pyfits.PrimaryHDU(data)
hdu.writeto(outfile, clobber=True)

