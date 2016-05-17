import Graffity
import os
import glob
import matplotlib.pyplot as pyplot
import numpy

fig = pyplot.figure(0)
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

datadir = '/home/deen/Data/GRAVITY/MATLAB_DATA/2016-03-15/ANAMORPH-193144/'
modalBasis = '/home/deen/Data/GRAVITY/MATLAB_DATA/2016-03-15/CDMS/RecnOptimiser.ModalBasis.fits'

loops = glob.glob(datadir+'CIAO_LOOP*.fits')
loops.sort()

wfs = Graffity.WFS()
wfs.expose()

anamorphose = []

for loop in loops:
    CB = Graffity.CircularBuffer(df = loop, ModalBasis=modalBasis)
    CB.extractModulation()
    CB.extractIM()
    anamorphose.append(CB)

print 'All Anamorphose extracted'
for i in range(60):
    centroidResponses = []
    for angle in anamorphose:
        centroidResponses.append(angle.IM[i])
    centroidResponses = numpy.array(centroidResponses)

    ax.clear()
    ax.imshow(centroidResponses)
    ax.figure.show()
    raw_input()

