import scipy
import numpy
import matplotlib.pyplot as pyplot
import pyfits

steps = numpy.linspace(-0.4, 0.4, 301)

subsample = 50

disturb = []

for step in steps:
    frame = numpy.zeros(2)
    frame[0] = step
    for i in range(subsample*2):
        disturb.append(frame)

disturb = numpy.array(disturb, dtype=numpy.float32)
outfile = pyfits.PrimaryHDU(data=disturb)
outfile.writeto("TTM_scan.fits", clobber=True)
