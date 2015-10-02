import Graffity
import numpy
import scipy
import matplotlib.pyplot as pyplot

wave = 632.8

ciao = Graffity.WFS(wavelength=1800.0)

zern = [0.0, 0.0, 0.0, 0.0, 0.0]
pupil = [0.0, 0.0]
actPoke = numpy.zeros(60, dtype=numpy.float32)
derotAngle = 0.00
clockingAngle = 0.0

# Take the flat-wavefront image
ciao.setupInstrument(zern, pupil, actPoke, derotAngle, clockingAngle)
ciao.expose()

f0 = pyplot.figure(0)
f1 = pyplot.figure(1)
f2 = pyplot.figure(2)
f3 = pyplot.figure(3)
f4 = pyplot.figure(4)

ax0 = f0.add_axes([0.1, 0.1, 0.8, 0.8])
ax1 = f1.add_axes([0.1, 0.1, 0.8, 0.8])
ax2 = f2.add_axes([0.1, 0.1, 0.8, 0.8])
ax3 = f3.add_axes([0.1, 0.1, 0.8, 0.8])
ax4 = f4.add_axes([0.1, 0.1, 0.8, 0.8])

wferror = numpy.linspace(-3*wave, 3*wave, num=13)

clockingAngle = 0.02
for rms in wferror:
    print rms
    zern = [rms, 0.0, 0.0, 0.0, 0.0]
    ciao.setupInstrument(zern, pupil, actPoke, derotAngle, clockingAngle)
    ciao.expose()
    zern = [0.0, rms, 0.0, 0.0, 0.0]
    ciao.setupInstrument(zern, pupil, actPoke, derotAngle, clockingAngle)
    ciao.expose()
    zern = [0.0, 0.0, rms, 0.0, 0.0]
    ciao.setupInstrument(zern, pupil, actPoke, derotAngle, clockingAngle)
    ciao.expose()
    zern = [0.0, 0.0, 0.0, rms, 0.0]
    ciao.setupInstrument(zern, pupil, actPoke, derotAngle, clockingAngle)
    ciao.expose()
    zern = [0.0, 0.0, 0.0, 0.0, rms]
    ciao.setupInstrument(zern, pupil, actPoke, derotAngle, clockingAngle)
    ciao.expose()

centroids = numpy.array(ciao.centroids)
flat = centroids[0]
tip = numpy.average(numpy.abs(centroids[[i*5+1 for i in range(len(wferror))]]-flat), axis=1).transpose()
tilt = numpy.average(numpy.abs(centroids[[i*5+2 for i in range(len(wferror))]]-flat), axis=1).transpose()
focus = numpy.average(numpy.abs(centroids[[i*5+3 for i in range(len(wferror))]]-flat), axis=1).transpose()
astig1 = numpy.average(numpy.abs(centroids[[i*5+4 for i in range(len(wferror))]]-flat), axis=1).transpose()
astig2 = numpy.average(numpy.abs(centroids[[i*5+5 for i in range(len(wferror))]]-flat), axis=1).transpose()


ax0.plot(wferror, tip[0], '-o')
ax0.plot(wferror, tip[1], '-o')
ax0.set_title('Tip')
ax1.plot(wferror, tilt[0], '-o')
ax1.plot(wferror, tilt[1], '-o')
ax1.set_title('Tilt')
ax2.plot(wferror, focus[0], '-o')
ax2.plot(wferror, focus[1], '-o')
ax2.set_title('Focus')
ax3.plot(wferror, astig1[0], '-o')
ax3.plot(wferror, astig1[1], '-o')
ax3.set_title('Oblique Astigmatism')
ax4.plot(wferror, astig2[0], '-o')
ax4.plot(wferror, astig2[1], '-o')
ax4.set_title('Vertical Astigmatism')

f0.show()
f1.show()
f2.show()
f3.show()
f4.show()

f0.savefig('tip.png')
f1.savefig('tilt.png')
f2.savefig('focus.png')
f3.savefig('ObliqAstig.png')
f4.savefig('VertAstig.png')
