import scipy
import numpy
import pyfits
import matplotlib.pyplot as pyplot
import Graffity
import glob

fig = pyplot.figure(0)
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
ax.clear()

datadir = '/home/deen/Data/GRAVITY/MATLAB_DATA/2016-06-16/'

datafiles = glob.glob(datadir+'DATA_LOGGER*')
#'DATA_LOGGER-080224/CIAO_LOOP_0001.fits'

voltages = []
enabled = []
high = []

for df in datafiles:
    loop = Graffity.LoopAnalyzer(df=df+'/CIAO_LOOP_0001.fits')
    loop.predict()
    high.append(loop.CB.header.get('HIERARCH ESO AOS HOCTR SMA HIGH'))
    sma_enabled = loop.CB.header.get('HIERARCH ESO AOS HOCTR SMA ENABLE')
    if sma_enabled:
        enabled.append(1)
    else:
        enabled.append(0)
    voltages.append(numpy.mean((loop.CB.HODM*loop.CB.HODM)**0.5)*400)

ax.plot(high, voltages)
#ax.set_aspect('auto')
fig.show()
