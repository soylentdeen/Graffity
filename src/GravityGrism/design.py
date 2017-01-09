import scipy
import numpy
from matplotlib import pyplot
import Grism

fig = pyplot.figure(0)
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
ax.clear()

G = Grism.Grism(n=3.43)
C = Grism.CCD()

focalLength = []
slitImage = []
slitWidth = []

"""
for f1 in numpy.arange(2.0, 50.0):
    for npix in numpy.arange(2.0, 100.0):
        focalLength.append(f1)
        slitImage.append(npix)
        Instrument = Grism.GrismInstrument(Grism=G, f1=f1)
        Instrument.optimize(npix=npix)
        slitWidth.append(Instrument.Grism.deltaX_slit)


focalLength = numpy.array(focalLength)
slitImage = numpy.array(slitImage)
slitWidth = numpy.array(slitWidth)
"""

Instrument = Grism.GrismInstrument(Grism=G, f1=54.6, f2=134.3, CCD=C)
Instrument.optimize(npix=2.0, Lambda = 2200.0, dLambda=500.0, order=5)

print("Sigma = %.3f, %.1f lines/mm" % (Instrument.Grism.sigma, 1.0/(Instrument.Grism.sigma/1000.0)))
print("Resolving Power = %.3f" % Instrument.ResolvingPower)
print("Delta = %.3f" % Instrument.Grism.delta)
print("Slit Size = %.3f" % Instrument.deltaX_slit)


#ax.scatter(focalLength, slitWidth, c = slitImage)
#fig.show()




