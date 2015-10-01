import Graffity
import numpy
import scipy
import matplotlib.pyplot as pyplot

ciao = Graffity.WFS()

zern = [0.0, 0.0, 0.0, 0.0, 0.0]
pupil = [0.0, 0.0]
actPoke = numpy.zeros(60, dtype=numpy.float32)
angle = 0.0

ciao.setupInstrument(zern, pupil, actPoke, angle)

fig = pyplot.figure(0)
fig.clear()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
ax.matshow(ciao.detector.z[0])
fig.show()
