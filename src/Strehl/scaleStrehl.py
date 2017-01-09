import scipy
from matplotlib import pyplot
import numpy

WFE = numpy.arange(0.0, 0.5, 0.05)

hband = 1.6
kband = 2.2

Strehl_Hband = numpy.exp(-(2.0*numpy.pi*WFE/hband)**2.0)
Strehl_Kband = numpy.exp(-(2.0*numpy.pi*WFE/kband)**2.0)

fig = pyplot.figure(0)
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
ax.clear()

ax.plot(WFE, Strehl_Hband)
ax.plot(WFE, Strehl_Kband)

fig.show()

