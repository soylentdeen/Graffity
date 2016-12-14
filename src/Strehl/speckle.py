import scipy
import numpy
from matplotlib import pyplot
import Graffity
import astropy.io.fits as pyfits
import glob
import time

fig = pyplot.figure(0)
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
ax.clear()

datadir = '/home/cdeen/Data/GRAVITY/2016-09-21/'

files = glob.glob(datadir+'*.cropped')
files.sort()

PSF = Graffity.PSF()
PSF.generateOTF()
theoreticalSum = numpy.sum(PSF.getPSF())

SRs = {}
SRs[0] = []
SRs[1] = []
SRs[2] = []
SRs[3] = []

dSRs = {}
dSRs[0] = []
dSRs[1] = []
dSRs[2] = []
dSRs[3] = []
timestamp = []

out = open('strehlRatios.dat', 'w')

for f in files[:16]:
    timestamp.append(time.mktime(time.strptime(f.split('GRAVITY.')[1].split('.fits.')[0], '%Y-%m-%dT%H-%M-%S')))
    AcqCamImage = Graffity.AcqCamImage(datadir=datadir, df=f.split('/')[-1])

    AcqCamImage.findPeaks(size=60)
    AcqCamImage.stack()

    StrehlRatios = {}
    for ut in range(4):
        StrehlRatios[ut] = []
        for frame, PSF in zip(AcqCamImage.speckleStamps[ut], AcqCamImage.postageStamps[ut]):
            
            for star in frame:
                ax.clear()
                ax.matshow(star - PSF)
                print numpy.max(star-PSF)
                fig.show()
                raw_input()

        StrehlRatios[ut] = numpy.array(StrehlRatios[ut])
        SRs[ut].append(numpy.mean(StrehlRatios[ut]))
        dSRs[ut].append(numpy.std(StrehlRatios[ut]))

    out.write("%d %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n" % (timestamp[-1], SRs[0][-1], dSRs[0][-1],
                  SRs[1][-1], dSRs[1][-1], SRs[2][-1], dSRs[2][-1], SRs[3][-1], dSRs[3][-1]))

out.close()

ax.scatter(timestamp, SRs[0], color = 'b')
ax.scatter(timestamp, SRs[1], color = 'g')
ax.scatter(timestamp, SRs[2], color = 'r')
ax.scatter(timestamp, SRs[3], color = 'y')
fig.show()

