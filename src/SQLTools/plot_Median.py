import numpy
import scipy
import matplotlib.pyplot as pyplot
import astropy.io.fits as pyfits

CIAO1 = pyfits.getdata('CIAO1.fits')
CIAO2 = pyfits.getdata('CIAO2.fits')
CIAO3 = pyfits.getdata('CIAO3.fits')
CIAO4 = pyfits.getdata('CIAO4.fits')

fig = pyplot.figure(0)
fig.clear()
ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.4])
ax2 = fig.add_axes([0.1, 0.5, 0.8, 0.4])

ax1.plot([numpy.log10(2.2), numpy.log10(2.2)], [-19, -11], color = 'm')
ax2.plot([numpy.log10(2.2), numpy.log10(2.2)], [-19, -11], color = 'm')
ax1.plot(CIAO1[0], CIAO1[1], color = 'b', label = 'CIAO1')
ax2.plot(CIAO1[0], CIAO1[2], color = 'b')
ax1.plot(CIAO2[0], CIAO2[1], color = 'g', label = 'CIAO2')
ax2.plot(CIAO2[0], CIAO2[2], color = 'g')
ax1.plot(CIAO3[0], CIAO3[1], color = 'r', label = 'CIAO3')
ax2.plot(CIAO3[0], CIAO3[2], color = 'r')
ax1.plot(CIAO4[0], CIAO4[1], color = 'k', label = 'CIAO4')
ax2.plot(CIAO4[0], CIAO4[2], color = 'k')

ax1.text(0.5, -18, "Tip")
ax2.text(0.5, -18, "Tilt")
ax1.legend()
fig.show()
fig.savefig("MedianTipTiltPower.png")
