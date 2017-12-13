import scipy
import numpy
import sys
import glob
sys.path.append('../')
import Graffity
from matplotlib import pyplot
import astropy.io.fits as pyfits
from scipy.optimize import leastsq



fig = pyplot.figure(1)
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
ax.clear()

fig2 = pyplot.figure(2)
axes = []
axes.append(fig2.add_axes([0.1, 0.1, 0.4, 0.4]))
axes.append(fig2.add_axes([0.1, 0.5, 0.4, 0.4]))
axes.append(fig2.add_axes([0.5, 0.1, 0.4, 0.4]))
axes.append(fig2.add_axes([0.5, 0.5, 0.4, 0.4]))

fig3 = pyplot.figure(3)
ax3 = fig3.add_axes([0.1, 0.1, 0.8, 0.8])
ax3.clear()

datadir = '/home/cdeen/Data/GRAVITY/AcqCam/'
df = glob.glob(datadir+'Calib*.fits')
#df = glob.glob(datadir+'GRAVI*2017-08-10T02:21:54.904_dualscip2vmred.fits')

subImageSize = 250

Telescopes= [1, 2, 3, 4]
colors = ["g", "b", "r", "y"]

Strehl = {}
Strehl[1] = []
Strehl[2] = []
Strehl[3] = []
Strehl[4] = []

def binStrehl(SR, AT, FT):
    new = []
    startTime = FT[0] - numpy.mean(numpy.diff(FT))
    stopTime = (FT[1] + FT[0])/2.0
    included = (startTime < AT) & (stopTime > AT)
    new.append(numpy.mean(SR[included]))
    for t in range(len(FT)-2):
        startTime = (FT[t] +FT[t+1])/2.0
        stopTime = (FT[t+1] + FT[t+2])/2.0
        included = (startTime < AT) & (stopTime > AT)
        new.append(numpy.mean(SR[included]))
    startTime = (FT[-2] - numpy.mean(numpy.diff(FT)))/2.0
    stopTime = FT[-1]
    included = (startTime < AT) & (stopTime > AT)
    new.append(numpy.mean(SR[included]))

    return numpy.array(new)

def newImg(p, x, y):
    #return p[0]*numpy.exp(-(x-p[1])**2.0/(numpy.cos(p[5])*p[2]+numpy.sin(p[5])*p[4]) - ((y-p[3])**2.0/(-numpy.sin(p[5])*p[2]+numpy.cos(p[5])*p[4])))
    return p[0]*numpy.exp(-(x-p[1])**2.0/(p[2]) - ((y-p[3])**2.0/(p[4])))

def calcEllipse(p, img, xi, yi):
    return (newImg(p, xi, yi).flatten() - img.flatten())**2.0
    
    
    
def fitEllipse(img):
    c = 1.0
    rx = 1.0
    ry = 1.0
    xc = 20.0
    yc = 20.0
    theta = 0.0
    x = numpy.arange(img.shape[0])
    y = numpy.arange(img.shape[1])
    xi, yi = numpy.meshgrid(x, y)
    fit = leastsq(calcEllipse, [c, xc, rx, yc, ry], args=(img, xi, yi))
    ratio = fit[0][2]/fit[0][4]
    return ratio

Perfect = Graffity.PSF(sizeInPix=41)
Perfect.generateOTF()

ImageStack = {}
ImageStack[1] = []
ImageStack[2] = []
ImageStack[3] = []
ImageStack[4] = []
elipse = {}
elipse[1] = []
elipse[2] = []
elipse[3] = []
elipse[4] = []
AcqFlux = {}
AcqFlux[1] = []
AcqFlux[2] = []
AcqFlux[3] = []
AcqFlux[4] = []

for f in df:
    AcqCamImages = pyfits.getdata(f)
    for im in AcqCamImages:
        for i, xcoord, ycoord in zip(Telescopes, [124, 367, 618, 873], [86, 86, 95, 90]):
            AcqFlux[i].append(numpy.max(im[ycoord-1:ycoord+1,xcoord-1:xcoord+1]))
            ImageStack[i].append(im[ycoord-20:ycoord+19, xcoord-20:xcoord+19])
            ImageStack[i][-1] /= numpy.max(ImageStack[i][-1])
            elipse[i].append(fitEllipse(ImageStack[i][-1]))
            strehl = Perfect.calcStrehl(im[ycoord-20:ycoord+19, xcoord-20:xcoord+19])
            if (numpy.isnan(strehl)):
                print asdf
            #print("Tel %d: X = %.2f, Y = %.2f, Strehl = %.3f" % (i, xcoord, ycoord, strehl))
            Strehl[i].append(strehl)

binnedStrehl = {}
binnedEllipse = {}
binnedAcqFlux = {}
xcross_section = {}
ycross_section = {}
largestFlux = 0
for i in Telescopes:
    #ImageStack[i] = Perfect.normalize(numpy.mean(numpy.array(ImageStack[i]), axis=0))
    #ImageStack[i] = ImageStack[i] * numpy.sum(Perfect.PSF) / numpy.sum(ImageStack[i])
    #ImageStack[i] = numpy.std(numpy.array(ImageStack[i]), axis=0)
    Strehl[i] = numpy.array(Strehl[i])
    elipse[i] = numpy.array(elipse[i])
    AcqFlux[i] = numpy.array(AcqFlux[i])
    #ax.scatter(binnedEllipse[i], FTFlux[i-1], c=colors[i-1])
    ax.scatter(Strehl[i], AcqFlux[i], c=colors[i-1])
    #ax.scatter(binnedStrehl[i], FTFlux[i-1], c=colors[i-1])
    #axes[i-1].matshow(ImageStack[i])
    axes[i-1].clear()
    xcross_section[i] = numpy.sum(ImageStack[i], axis=1)
    if numpy.sum(xcross_section[i]) > largestFlux:
        largestFlux = numpy.sum(xcross_section[i])
    ycross_section[i] = numpy.sum(ImageStack[i], axis=0)
    if numpy.sum(ycross_section[i]) > largestFlux:
        largestFlux = numpy.sum(ycross_section[i])
    axes[i-1].matshow(numpy.median(ImageStack[i], axis=0))
    #axes[i-1].plot(xcross_section[i], c=colors[i-1])
    #axes[i-1].plot(ycross_section[i], c=colors[i-1])
    #axes[i-1].plot(numpy.sum(Perfect.PSF, axis=1), c='k')



for i in Telescopes:
    ax3.plot(xcross_section[i], c=colors[i-1])
    ax3.plot(ycross_section[i], c=colors[i-1])

ax3.plot(numpy.sum(Perfect.PSF, axis=1), c='k')

ax.set_xlabel("H-band Strehl")
ax.set_ylabel("Peak Flux")
ax.set_title("Strehl on Calibration Images")
fig.show()
fig.savefig("Calibration_Strehl.png")
fig2.show()
fig3.show()
