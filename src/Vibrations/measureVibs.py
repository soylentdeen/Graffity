import Graffity
from matplotlib import pyplot
import scipy
import numpy

fig = pyplot.figure(0)
ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.2])
ax1.clear()
ax2 = fig.add_axes([0.1, 0.3, 0.8, 0.2])
ax2.clear()
ax3 = fig.add_axes([0.1, 0.5, 0.8, 0.2])
ax3.clear()
ax4 = fig.add_axes([0.1, 0.7, 0.8, 0.2])
ax4.clear()

fig2 = pyplot.figure(1)
ax = fig2.add_axes([0.1, 0.1, 0.8, 0.8])
ax.clear()


datadir = 'CIAO_1/2017-08/2017-08-09/DATA_LOGGER-030523'

dl = Graffity.DataLogger(datadir)
dl.loadData()
dl.computeStrehl()
dl.measureVibs(frequencies=[124.0])


scale = 1e9*numpy.sqrt(dl.ZPowerdFrequencies)
for i in [0,1,2,7,8]:
    ax4.plot(dl.ZPowerFrequencies, dl.ZPowerSlopes[i]*scale)
    ax3.plot(dl.ZPowerFrequencies, numpy.cumsum(dl.ZPowerSlopes[i]*scale))
    ax2.plot(dl.ZPowerFrequencies, numpy.sqrt(dl.ZPowerCommands[i]*scale))
    ax1.plot(dl.ZPowerFrequencies, numpy.cumsum(numpy.sqrt(dl.ZPowerCommands[i]*scale)))



fig.show()

for key in dl.vibPower.keys():
    ax.scatter(dl.vibPower[key]['Freq'], numpy.log10(dl.vibPower[key]['CommPower']))

fig2.show()
