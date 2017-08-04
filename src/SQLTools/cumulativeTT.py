import astropy.io.fits as pyfits
import matplotlib.pyplot as pyplot
import numpy
import scipy

CIAO1 = pyfits.getdata("CIAO1.fits")
CIAO2 = pyfits.getdata("CIAO2.fits")
CIAO3 = pyfits.getdata("CIAO3.fits")
CIAO4 = pyfits.getdata("CIAO4.fits")

fig = pyplot.figure(0)
ax1 = fig.add_axes([0.1, 0.1, 0.4, 0.4])
ax2 = fig.add_axes([0.1, 0.5, 0.4, 0.4])
ax3 = fig.add_axes([0.5, 0.1, 0.4, 0.4])
ax4 = fig.add_axes([0.5, 0.5, 0.4, 0.4])

for ciao, ax in zip([CIAO1, CIAO2, CIAO3, CIAO4], [ax1, ax2, ax3, ax4]):
    freq = ciao[0]
    tip = ciao[1]
    tipCDF = [tip[0]]
    for t in tip[1:]:
        tipCDF.append(t+tipCDF[-1])
    tilt = ciao[2]
    tiltCDF = [tilt[0]]
    for t in tilt[1:]:
        tiltCDF.append(t+tiltCDF[-1])

    ax.step(freq, tipCDF, color = 'b')
    ax.step(freq, tiltCDF, color = 'r')


fig.show()

