import Graffity
import numpy
import matplotlib.pyplot as pyplot
import scipy
import astropy.io.fits as pyfits
import sys
import glob
from scipy import optimize

def fitSlope(tt, centroid):
    errfunc_X = lambda p, tt : p[0] + p[2]*tt[:,0]+p[4]*tt[:,1]
    errfunc_Y = lambda p, tt : p[1] + p[3]*tt[:,0]+p[5]*tt[:,1]
    #fitfunc = lambda p, x, y, th : numpy.abs(errfunc_X(p, th) - x) + numpy.abs(errfunc_Y(p, th) - y)
    fitfunc = lambda p, x, y, tt : numpy.abs(errfunc_X(p, tt) - x) + numpy.abs(errfunc_Y(p, tt) - y)
    x = centroid[:,0]
    y = centroid[:,1]
    params = [3.5, 3.5, 0.5, 0.05, 0.5, 0.15]
    pfit = optimize.leastsq(fitfunc, params, args = (x, y, tt), full_output=True)
    return pfit

def extractSlope(tt, centroids):
    matrix = scipy.linalg.pinv(tt).dot(centroids)
    return matrix

datadir = '/Users/ciao/Data/CIAO/DATA/2016-05-17_4/FOV-153439/'

loops = glob.glob(datadir+'CIAO_LOOP*.fits')
loops.sort()
pixels = glob.glob(datadir+'CIAO_PIXELS*.fits')
pixels.sort()


pyplot.rcParams['font.size'] = 18
fig = pyplot.figure(0)
fig.clear()
ax = fig.add_axes([0.15, 0.15, 0.8, 0.8])
fig2 = pyplot.figure(2)
fig2.clear()
ax2 = fig2.add_axes([0.15, 0.15, 0.8, 0.8])

tip = []
tilt = []
xcentroids = []
ycentroids = []

npts = len(loops)-1
delta = 2.0*numpy.pi/npts
theta = []

WFS = Graffity.WFS()


for pixel, loop in zip(pixels, loops):

    loop_data = pyfits.getdata(loop)
    pixel_data = pyfits.getdata(pixel)
    
    TTM = loop_data.field("ITTM_Positions")
    loop_frames = loop_data.field('FrameCounter')
    pixel_frames = pixel_data.field('FrameCounter')
    images = pixel_data.field('Pixels')
    
    centroids = []
    tiptilt = []

    ax.clear()
    ax2.clear()

    ft = numpy.fft.fftshift(numpy.fft.fft2(images[0].reshape(72,72)))
    ax.imshow(ft.real)
    ax2.imshow(ft.imag)

    fig.show()
    fig2.show()

    blah = input("Enter to continue")

    ax.clear()
    ax2.clear()

    ft = numpy.fft.fftshift(numpy.fft.fft2(images[150].reshape(72,72)))
    ax.imshow(ft.real)
    ax2.imshow(ft.imag)

    fig.show()
    fig2.show()

    blah = input("Enter to continue")


    for im, f in zip(images, pixel_frames):
        ax.clear()
        ft = numpy.fft.fftshift(numpy.fft.fft2(im.reshape(72,72)))
        ax.imshow(ft.imag)
        fig.show()
        blah = input()
        """
        try:
            tiptilt.append(TTM[loop_frames == f][0])
            centroids.append(WFS.detector.measureCentroids(im))

        except:
            pass
        #"""
    
    centroids = numpy.array(centroids)
    tiptilt = numpy.array(tiptilt)
    tiptilt -= tiptilt[0]
    tiptilt *= 5.0
    spots = numpy.average(centroids, axis=1)
    spots -= spots[0]
    PlateScale = extractSlope(tiptilt, spots)
    print(PlateScale)

    ax.plot(spots[:,0], spots[:,1])
    values = tiptilt.dot(PlateScale)
    ax.scatter(values[:,0], values[:,1])
    fig.show()
    input()

fig.show()


centx = circle[0]
centy = circle[1]
rx = circle[2]
ry = circle[3]
points_x = []
points_y = []

newTheta = numpy.linspace(0, numpy.pi*2.0)

for th in newTheta:
    points_x.append(centx + rx*numpy.cos(th))
    points_y.append(centy + ry*numpy.sin(th))


ax.scatter(tip, tilt, s=30)
ax.plot(points_x, points_y, color = 'r', lw=3.0)

plateScale = numpy.mean([numpy.abs(rx), numpy.abs(ry)])

ax.text(0.2, 0.9, 'Plate Scale: %.3f TTMU/pix' % plateScale, 
        fontsize=20, transform=ax.transAxes)

ax.set_ylabel("Tilt Position")
ax.set_xlabel("Tip Position")

ax.set_aspect('equal')
#ax.scatter(x, y)
#ax.set_xbound(-0.1, 0.1)
#ax.set_ybound(-0.1, 0.1)
fig.show()
fig.savefig("PlateScale.png")
