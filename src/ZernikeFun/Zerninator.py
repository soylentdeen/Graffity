import scipy
import numpy
import matplotlib.pyplot as pyplot
import pyfits
import Graffity


class Map( object ):
    def __init__(self, default, ZernikeMatrix):
        self.default = default
        self.ZernikeMatrix = ZernikeMatrix
        self.result = default
        self.WFS = Graffity.WFS()
        #self.LA = Graffity.lensletArray()

    def applyCoefficients(self, coefficients):
        self.WFS.setupInstrument(zern=coefficients)
        self.WFS.expose()

    def display(self, ax):
        ax.clear()
        #for i in range(len(self.LA.coordinates)):
        #    nominalX = self.LA.coordinates[i][0]/24.0
        #    nominalY = self.LA.coordinates[i][1]/24.0
        #    ax.arrow(nominalX, nominalY, self.result[i*2], 
        #            self.result[i*2+1], head_width=0.5, head_length=1.0)
        ax.figure.show()
        print numpy.max(numpy.abs(self.result))
        raw_input()

    def save(self, name):
        hdu = pyfits.PrimaryHDU(data = self.result)
        hdu.writeto(name, clobber=True)

fig = pyplot.figure(0)
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
ax.clear()

#defaultRefSlopes = pyfits.getdata('../../data/Acq.DET1.REFSLP.fits')
defaultRefSlopes = numpy.ones(136)*3.5
defaultActPos = pyfits.getdata('../../data/HOCtr.ACT_POS_REF_MAP.fits')

#S2Z = pyfits.getdata('../../data/Slopes2Z.fits')
#Z2S = numpy.linalg.pinv(S2Z)/400.0
#Z2S = pyfits.getdata('../../data/Z2S_136s_119z.fits').T
#Z2S = pyfits.getdata('../../data/Slopes2Zer_CIAO.fits')
Z2S = pyfits.getdata('../../data/Zer2Slopes_CIAO.fits').T
Z2DM = pyfits.getdata('../../data/cimdatZernike2DM.fits')

RefSlopes = Map(defaultRefSlopes, Z2S)
DMPositions = Map(defaultActPos, Z2DM)

for i in range(12):
    for rms in numpy.linspace(-100, 100, 9):
        zernikes = numpy.zeros(20)
        if i < 2:
            zernikes[i] = rms/50.0
        else:
            zernikes[i] = rms/100.0
        RefSlopes.applyCoefficients(zernikes)
        name = './RefSlopes/RefSlopes_Z%2d_RMS%4d.fits' % (i, rms)
        RefSlopes.display(ax)
        RefSlopes.save(name)

        


#fig.show()
