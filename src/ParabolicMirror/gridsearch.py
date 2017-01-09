import scipy
import numpy
import matplotlib.pyplot as pyplot
import astropy.io.fits as pyfits
import glob


datadir = '/home/cdeen/Data/CIAO/UT3/Alignment/2016-08-31_3/PARABOLA_SEARCH-190246/'

files = glob.glob(datadir+'*.fits')

scan = []
flux = []
scans = []

for f in files:
    data = pyfits.getdata(f).field('Intensities')
    header = pyfits.getheader(f)
    scan.append(header.get('ESO TPL EXPNO'))
    flux.append(numpy.max(numpy.mean(data, axis=1)))
    scans.append(numpy.mean(data, axis=1))

middle = numpy.array([1717503, 1667465])

corner = middle - numpy.array([250000, 250000])

step = 5000

scan = numpy.array(scan)
order = numpy.argsort(scan)
scan = scan[order]*step+corner[0]
scans = numpy.array(scans)[order]
flux = numpy.array(flux)[order]

fig = pyplot.figure(0)
fig.clear()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

for s in scans[40:50]:
    ax.plot(s)
#ax.plot(scan, flux)

ax.set_ylabel("Average Subaperture Intensity")
ax.set_xlabel("Frame number")
ax.set_xbound(3200, 3500)

fig.show()
fig.savefig("Scan40to50.png")

raw_input()
ax.clear()
ax.plot(scan, flux)
ax.set_xlabel('Scan number')
ax.set_ylabel('Brightest Frame')
ax.get_xaxis().get_major_formatter().set_useOffset(False)
fig.show()
