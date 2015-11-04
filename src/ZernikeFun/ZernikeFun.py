import scipy
import numpy
#import matplotlib.pyplot as pyplot
import pyfits
import py_compile
import sys

#fig = pyplot.figure(0)
#fig.clear()
#ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

Z2DM = pyfits.getdata('cimdatZernike2DM.fits')

nominal = pyfits.getdata('QSPATTERN.fits')

outName = 'result.fits'

input_coeffs = numpy.array(sys.argv)

coeffs = numpy.zeros(len(Z2DM))

coeffs[0:len(input_coeffs)-1] = numpy.array(input_coeffs[1:])

delta = numpy.dot(coeffs*1e-6, Z2DM)

output = nominal + delta

actnum = numpy.arange(60)

out = pyfits.PrimaryHDU(output)
out.writeto(outName, clobber=True)

#ax.plot(actnum, nominal[0])
#ax.plot(actnum, output[0])
#fig.show()
