import Graffity
import scipy
import numpy
import matplotlib.pyplot as pyplot

wave = 632.8

ciao = Graffity.WFS(wavelength = 1800.0)

clockingAngles = numpy.linspace(-0.10, 0.10, num=21)

zern = [0.0, 0.0, 0.0, 0.0, 0.0]
pupil = [0.0, 0.0]
actPoke = numpy.zeros(60, dtype=numpy.float32)
derotAngle = 0.00
clockingAngle = 0.0

# Take the flat-wavefront image
ciao.setupInstrument(zern, pupil, actPoke, derotAngle, clockingAngle)
ciao.expose()

for cA in clockingAngles:
    print(cA)
    ciao.setupInstrument(zern, pupil, actPoke, derotAngle, cA)
    ciao.expose()

centroids = numpy.array(ciao.centroids)
flat = centroids[0]

slopes = numpy.max(centroids[1:]-flat, axis=1)

fig = pyplot.figure(0)
fig.clear()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
ax.plot(clockingAngles, slopes[:,0], color='b', marker = 'o')
ax.plot(clockingAngles, slopes[:,1], color='r', marker = 'o')
ax.set_xlabel("Clocking Angle (Degrees)")
ax.set_ylabel("Maximum Slope (Pixels)")
ax.set_title("Clocking angle dependence on Static Slopes")
fig.show()
fig.savefig("clocking.png")
