import scipy
import numpy
import Graffity
import matplotlib.pyplot as pyplot

fig = pyplot.figure(0)
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

#Number of rows, columns
nx = 21
ny = 20

datadir = '/home/deen/Data/GRAVITY/Alignment/29Oct2015/'

dark = Graffity.NGCImage(datadir+'dark.fits')
light = Graffity.NGCImage(datadir+'light.fits')

light.subtract_background(dark)
light.findSubapertureCenters(nx=nx, ny=ny)

light.findCentroids()

xangles, yangles = light.findAngles(ax=ax)

fig.show()
