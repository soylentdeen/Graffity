import scipy
import numpy
from matplotlib import pyplot
import Graffity
import astropy.io.fits as pyfits

fig = pyplot.figure(0)
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
ax.clear()

AcqCamImage = Graffity.AcqCamImage(df='/data/cdeen/Data/GRAVITY/ut2.ssa')

AcqCamImage.findPeaks(size=60)
AcqCamImage.stack()
PSF = Graffity.PSF()
PSF.generateOTF()

theoreticalSum = numpy.sum(PSF.getPSF())
measuredSum = numpy.sum(AcqCamImage.postageStamp)

strehlRatio = theoreticalSum/measuredSum
print strehlRatio

#ax.matshow(peaks)
ax.matshow(AcqCamImage.postageStamp)
#ax.matshow(PSF.getPSF())
fig.show()
