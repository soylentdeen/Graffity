import scipy
import numpy
import matplotlib.pyplot as pyplot
import pyfits

steps = numpy.linspace(-0.5, 0.5, 11)

disturb = []

for step in steps:
    frame = numpy.zeros(136)
    frame[::2] = step
    for i in range(2000):
        disturb.append(frame)

outfile = pyfits.PrimaryHDU(data=numpy.array(disturb, dtype=numpy.float32))
outfile.writeto("refSlopeDisturb.fits", clobber=True)
