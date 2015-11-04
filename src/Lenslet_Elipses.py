import scipy
import numpy
import pyfits
import Graffity
import matplotlib.pyplot as pyplot



datdir= '/home/deen/Data/GRAVITY/Alignment/29Oct2015/'

dark = Graffity.NGCImage(datdir+'dark.fits')
light = Graffity.NGCImage(datdir+'light.fits')

light.subtract_background(dark)

light.findSubapertureCenters(21, 20)

cutout = light.extractCutout(10, 10)

fig1 = pyplot.figure(1)
fig1.clear()
ax1 = fig1.add_axes([0.1, 0.1, 0.8, 0.8])

fig2 = pyplot.figure(2)
fig2.clear()
ax2 = fig2.add_axes([0.1, 0.1, 0.8, 0.8])

fig3 = pyplot.figure(3)
fig3.clear()
ax3 = fig3.add_axes([0.1, 0.1, 0.8, 0.8])
"""
sig_x = 3.0
sig_y = 3.0

theta = numpy.deg2rad(45)

xsm = numpy.linspace(0, 8)
ysm = numpy.linspace(0, 8)

xcenter = 4.0
ycenter = 4.0

image = numpy.zeros((len(xsm), len(ysm)))

a = numpy.cos(theta)**2/(2.0*sig_x**2.0) + numpy.sin(theta)**2./(2.0*sig_y**2.0)
b = -numpy.sin(2.0*theta)**2/(4.0*sig_x**2.0) + numpy.sin(2.0*theta)**2./(4.0*sig_y**2.0)
c = numpy.sin(theta)**2/(2.0*sig_x**2.0) + numpy.cos(theta)**2./(2.0*sig_y**2.0)

Amplitude = 1.0

for x in range(len(xsm)):
    for y in range(len(ysm)):
        image[x][y] = Amplitude*numpy.exp(-(a*(xsm[x]-xcenter)**2.0 +
            2.0*b*(xsm[x]-xcenter)*(ysm[y]-ycenter) + c*(ysm[y]-ycenter)**2.0))

pixels = numpy.zeros((8, 8))
for x in range(len(pixels)):
    for y in range(len(pixels[0])):
        goodx = numpy.arange(len(xsm))[(x <= xsm) & (xsm < x+1)]
        goody = numpy.arange(len(ysm))[(y <= ysm) & (ysm < y+1)]
        pixels[x][y] = numpy.sum(image[goodx][:, goody])

for x in range(len(xsm)):
    for y in range(len(ysm)):
        image[x][y] = 1.0/((1.0 +numpy.abs((xsm[x]-xcenter)/sig_x)**5.0) * 
                (1.0+numpy.abs((ysm[y]-ycenter)/sig_y)**5.0))

ax2.matshow(image)
#"""
#ax.matshow(light.subtracted)
centers = numpy.zeros((light.nx*light.ny, 2))
for x in range(light.nx):
    for y in range(light.ny):
        centers[x*light.ny+y, 0] = light.xcenters[x]
        centers[x*light.ny+y, 1] = light.ycenters[y]

#ax.scatter(centers[:,0], centers[:,1])
#ax.matshow(cutout)

#fit=light.fitElipse(cutout.ravel())
light.twirl()
#light.fitElipses(ax1=ax1, ax2=ax2, ax3=ax3)
#"""

ax1.matshow(light.subtracted)
ax2.matshow(light.subtracted)
junk = ax2.scatter(light.xc, light.yc, c=light.sig_x/light.sig_y, s= 50.0, vmin=1.0, vmax=2.0)

fig2.colorbar(junk)

fig1.show()
fig2.show()
