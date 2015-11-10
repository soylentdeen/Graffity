import scipy
import numpy
import Graffity

#Number of rows, columns
nx = 21
ny = 20

datadir = '/home/deen/Data/GRAVITY/Alignment/29Oct2015/'

dark = Graffity.NGCImage(datadir+'dark.fits')
light = Graffity.NGCImage(datadir+'light.fits')

light.subtract_background(dark)
light.findSubapertureCenters(nx=nx, ny=ny)

light.findCentroids()

light.findAngles()
