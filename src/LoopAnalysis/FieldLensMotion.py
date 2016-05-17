import Graffity
import numpy
import matplotlib.pyplot as pyplot
import scipy
import pyfits
import sys
import glob

datadirs = glob.glob('/home/deen/Data/GRAVITY/MATLAB_DATA/2015-12-21/DATA_LOGGER*')

fX = []
fY = []

sX = []
sY = []

fig = pyplot.figure(0)
fig.clear()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

datadirs.sort()

for datadir in datadirs:
    df = datadir+'/CIAO_LOOP_0001.fits'
    data = pyfits.getdata(df)
    header = pyfits.getheader(df)

    slopes = numpy.average(data.field(4), axis=0)

    fX.append(header.get('HIERARCH ESO INS FLDL X'))
    fY.append(header.get('HIERARCH ESO INS FLDL Y'))

    sX.append(numpy.median(slopes[::2]))
    sY.append(numpy.median(slopes[1::2]))
    print fX[-1], fY[-1], sX[-1], sY[-1]

ax.scatter(sX, sY)
ax.set_xbound(-0.1, 0.1)
ax.set_ybound(-0.1, 0.1)
fig.show()
fig.savefig("FieldLens_XYSpots.png")
