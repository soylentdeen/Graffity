import scipy
import numpy
import Graffity
import matplotlib.pyplot as pyplot

fig = pyplot.figure(0)
fig.clear()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

#Number of rows, columns
nx = 20
ny = 20

datadir='/home/deen/Data/GRAVITY/DetectorTests/MOVPE1/Iteration4/'
light = Graffity.NGCImage(datadir+'0001_Light-on_Double_DIT0.05_Bias4.5.fits')
dark = Graffity.NGCImage(datadir+'0001_Light-off_Double_DIT0.05_Bias4.5.fits')

#datadir = '/home/deen/Data/GRAVITY/Alignment/29Oct2015/'
#dark = Graffity.NGCImage(datadir+'dark.fits')
#light = Graffity.NGCImage(datadir+'light.fits')

light.subtract_background(dark)
light.findSubapertureCenters(nx=nx, ny=ny)

light.findCentroids()

#ax.matshow(light.subtracted
light.findAngles(ax=ax)

print("Angle of Rows: %.3f +/- %.3f degrees" % (numpy.mean(light.xangle), 
    numpy.std(light.xangle)))
print("Angle of Colunns: %.3f +/- %.3f degrees" % (numpy.mean(light.yangle),
    numpy.std(light.yangle)))

print("X offset = %.4f" % numpy.mean(light.residual_x))
print("Y offset = %.4f" % numpy.mean(light.residual_y))
#ax.set_aspect('auto')
fig.show()
