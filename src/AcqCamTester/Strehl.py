import astropy.io.fits as pyfits
from matplotlib import pyplot
import numpy
import scipy
import glob

class BCI_Data(object):
    def __init__(self, data=None, error=None):
        self.data = {}
        self.error = None
        if data != None:
            for i in range(4):
                self.data[i] = data[i::4]
                if error != None:
                    if self.error == None:
                        self.error = {}
                    self.error[i] = error[i::4]

    def __pow__(self, power):
        for i in self.data.keys():
            if self.error != None:
                self.error[i] = power * self.data[i]*self.error[i]
            self.data[i] = self.data[i]**power
        return self

    def __add__(self, other):
        for i in self.data.keys():
            if self.error != None:
                self.error[i] = numpy.sqrt(self.error[i]**2.0 +other.error[i]**2.0)
            self.data[i] = self.data[i] + other.data[i]
        return self

    def getMedian(self):
        retval = {}
        for i in self.data.keys():
            retval[i] = numpy.median(self.data[i])

        self.median = retval
        return retval

    def rebin(self, time, otherTime, debug=False):
        new = BCI_Data()
        if self.error != None:
            new.error = {}
        for i in self.data.keys():
            new.data[i] = []
            if self.error != None:
                new.error[i] = []
            startTime = otherTime.data[i][0] - numpy.mean(numpy.diff(otherTime.data[i]))
            stopTime = (otherTime.data[i][1] + otherTime.data[i][0]) / 2.0
            included = (startTime < time.data[i]) & (stopTime > time.data[i])
            if self.error != None:
                new.data[i].append(numpy.sum(self.data[i][included]/self.error[i][included])/numpy.sum(
                              1.0/self.error[i][included]))
                new.error[i].append(numpy.sum(self.error[i][included]**-2.0)**-0.5)
            else:
                new.data[i].append(numpy.mean(self.data[i][included]))
            for t in range(len(otherTime.data[i])-2):
                startTime = (otherTime.data[i][t]+otherTime.data[i][t+1])/2.0
                stopTime = (otherTime.data[i][t+1]+otherTime.data[i][t+2])/2.0
                included = (startTime < time.data[i]) & (stopTime > time.data[i])
                if self.error != None:
                    new.data[i].append(numpy.sum(self.data[i][included]/
                                       self.error[i][included])/numpy.sum(1.0/self.error[i][included]))
                    new.error[i].append(numpy.sum(self.error[i][included]**-2.0)**-0.5)
                    if debug:
                        print new.data[i][-1], new.error[i][-1], numpy.mean(self.data[i][included])
                        raw_input()
                else:
                    new.data[i].append(numpy.mean(self.data[i][included]))
            startTime = otherTime.data[i][-2] + numpy.mean(numpy.diff(otherTime.data[i]))/2.0
            stopTime = otherTime.data[i][-1]
            included = (startTime < time.data[i]) & (stopTime > time.data[i])
            if self.error != None:
                new.data[i].append(numpy.sum(self.data[i][included]/
                                   self.error[i][included])/numpy.sum(1.0/self.error[i][included]))
                new.error[i].append(numpy.sum(self.error[i][included]**-2.0)**-0.5)
            else:
                new.data[i].append(numpy.mean(self.data[i][included]))

            new.data[i] = numpy.array(new.data[i])
            if self.error != None:
                new.error[i] = numpy.array(new.error[i])
        return new

class AcqCamData( object ):
    def __init__(self, fiberData=None, fluxData=None):
        self.Time = BCI_Data(data=fiberData.field('TIME'))
        self.Strehl = BCI_Data(data=fiberData.field('FIELD_STREHL'))
        self.SC_X = BCI_Data(data=fiberData.field('FIELD_SC_X'), error=fiberData.field('FIELD_SC_XERR'))
        self.SC_Y = BCI_Data(data=fiberData.field('FIELD_SC_Y'), error=fiberData.field('FIELD_SC_YERR'))
        #self.SC_XERR = BCI_Data(fiberData.field('FIELD_SC_XERR'))
        #self.SC_YERR = BCI_Data(fiberData.field('FIELD_SC_YERR'))
        self.FT_X = BCI_Data(data=fiberData.field('FIELD_FT_X'), error=fiberData.field('FIELD_FT_XERR'))
        self.FT_Y = BCI_Data(data=fiberData.field('FIELD_FT_Y'), error=fiberData.field('FIELD_FT_YERR'))
        #self.FT_XERR = BCI_Data(fiberData.field('FIELD_FT_XERR'))
        #self.FT_YERR = BCI_Data(fiberData.field('FIELD_FT_YERR'))
        self.SCALE = BCI_Data(data=fiberData.field('FIELD_SCALE'), error=fiberData.field('FIELD_SCALEERR'))
        #self.SCALEERR = BCI_Data(fiberData.field('FIELD_SCALEERR'))
        self.SC_FIBER_DX = BCI_Data(data=fiberData.field('FIELD_SC_FIBER_DX'), error=fiberData.field('FIELD_SC_FIBER_DXERR'))
        #self.SC_FIBER_DXERR = BCI_Data(fiberData.field('FIELD_SC_FIBER_DXERR'))
        self.SC_FIBER_DY = BCI_Data(data=fiberData.field('FIELD_SC_FIBER_DY'), error=fiberData.field('FIELD_SC_FIBER_DYERR'))
        #self.SC_FIBER_DYERR = BCI_Data(fiberData.field('FIELD_SC_FIBER_DYERR'))
        #self.FIBER_DELTA = (self.SC_FIBER_DX**2.0 + self.SC_FIBER_DY**2.0)**0.5
        #self.FIBER_DELTA_ERR = (self.SC_FIBER_DXERR**2.0 + self.SC_FIBER_DYERR**2.0)**0.5

        self.BCI_Time = BCI_Data(fluxData.field('TIME'))
        self.TOTALFLUX_SC = BCI_Data(fluxData.field('TOTALFLUX_SC'))
        self.TOTALFLUX_FT = BCI_Data(fluxData.field('TOTALFLUX_FT'))

    def binData(self):
        self.newStrehl = self.Strehl.rebin(self.Time, self.BCI_Time)
        self.newSC_X = self.SC_X.rebin(self.Time, self.BCI_Time)
        self.newSC_Y = self.SC_Y.rebin(self.Time, self.BCI_Time)
        self.newFT_X = self.FT_X.rebin(self.Time, self.BCI_Time)
        self.newFT_Y = self.FT_Y.rebin(self.Time, self.BCI_Time)
        self.newSCALE = self.SCALE.rebin(self.Time, self.BCI_Time)
        self.newSC_FIBER_DX = self.SC_FIBER_DX.rebin(self.Time, self.BCI_Time)
        self.newSC_FIBER_DY = self.SC_FIBER_DY.rebin(self.Time, self.BCI_Time)
        

    def findCorrelations(self, ax=None):
        colors = ['b', 'g', 'r', 'c']
        ax.clear()
        linex = numpy.array([numpy.ones(50), numpy.linspace(0, 0.15)])
        for i in range(4):
            A = []
            x = []
            for SR, flux in zip(self.newStrehl.data[i], self.TOTALFLUX_FT.data[i]):
                if not(numpy.isnan(SR)):
                    A.append([1.0, SR])
                    x.append(flux / self.TOTALFLUX_FT.median[i])
                    if ax != None:
                        ax.scatter([SR], [flux/self.TOTALFLUX_FT.median[i]], c=colors[i])

            A = numpy.array(A)
            x = numpy.array(x)
            B = numpy.linalg.pinv(A)
            fit = B.dot(x)
            if ax != None:
                ax.plot(linex[1], fit.dot(linex), c=colors[i])
            print fit

    def calcMedian(self):
        medianTOTALFLUX_SC = self.TOTALFLUX_SC.getMedian()
        medianTOTALFLUX_FT = self.TOTALFLUX_FT.getMedian()
        medianSCALE = self.SCALE.getMedian()
        medianSC_FIBER_DX = self.SC_FIBER_DX.getMedian()
        medianSC_FIBER_DY = self.SC_FIBER_DY.getMedian()
        medianStrehl = self.Strehl.getMedian()

    def plot(self, axes=None):
        colors = ['b', 'g', 'r', 'c']
        medianTOTALFLUX_SC = self.TOTALFLUX_SC.getMedian()
        medianTOTALFLUX_FT = self.TOTALFLUX_FT.getMedian()
        medianSCALE = self.SCALE.getMedian()
        medianSC_FIBER_DX = self.SC_FIBER_DX.getMedian()
        medianSC_FIBER_DY = self.SC_FIBER_DY.getMedian()
        medianStrehl = self.Strehl.getMedian()
        for i in range(4):
            color = colors[i]
            for ax, name in zip(axes, ['SC', 'FLUX', 'SCALE', 'FIBER', 'STREHL']):
                if name == 'SC':
                    ax.errorbar(self.SC_X.data[i], self.SC_Y.data[i], xerr=self.SC_X.error[i], yerr=self.SC_Y.error[i], c=color)
                    ax.errorbar(self.FT_X.data[i], self.FT_Y.data[i], xerr=self.FT_X.error[i], yerr=self.FT_Y.error[i], c=color)
                elif name == 'FLUX':
                    ax.plot(self.BCI_Time.data[i], self.TOTALFLUX_SC.data[i]/medianTOTALFLUX_SC[i], c=color)
                    ax.plot(self.BCI_Time.data[i], self.TOTALFLUX_FT.data[i]/medianTOTALFLUX_FT[i], c=color)
                elif name == 'SCALE':
                    #ax.errorbar(self.Time.data[i], self.SCALE.data[i], yerr=self.SCALE.error[i], c=color)
                    ax.errorbar(self.BCI_Time.data[i], self.newSCALE.data[i]-medianSCALE[i], yerr=self.newSCALE.error[i], c=color)
                elif name == 'FIBER':
                    #ax.errorbar(self.Time.data[i], self.FIBER_DELTA.data[i], yerr=self.FIBER_DELTA_ERR.data[i], c=color)
                    #ax.errorbar(self.Time.data[i], self.SC_FIBER_DX.data[i], yerr=self.SC_FIBER_DX.error[i], c=color, ls='-')
                    #ax.errorbar(self.Time.data[i], self.SC_FIBER_DY.data[i], yerr=self.SC_FIBER_DY.error[i], c=color, ls='-.')
                    ax.errorbar(self.BCI_Time.data[i], self.newSC_FIBER_DX.data[i], yerr=self.newSC_FIBER_DX.error[i], c=color, ls='-')
                    ax.errorbar(self.BCI_Time.data[i], self.newSC_FIBER_DY.data[i], yerr=self.newSC_FIBER_DY.error[i], c=color, ls='-.')
                    #ax.plot(self.Time.data[i], self.SC_FIBER_DX.data[i], c=color)
                    #ax.plot(self.Time.data[i], self.SC_FIBER_DY.data[i], c=color)
                    #ax.plot(self.BCI_Time.data[i], self.newSC_FIBER_DX.data[i]-medianSC_FIBER_DX[i], c=color)
                    #ax.plot(self.BCI_Time.data[i], self.newSC_FIBER_DY.data[i]-medianSC_FIBER_DY[i], c=color)
                elif name == 'STREHL':
                    #ax.plot(self.Time.data[i], self.Strehl.data[i], c=color)
                    ax.plot(self.BCI_Time.data[i], self.newStrehl.data[i]/medianStrehl[i], c=color)


df = '/home/cdeen/Data/GRAVITY/2017-08/2017-08-09/reduced/GRAVI.2017-08-10T00:35:33.638_dualcalp2vmred.fits'


def findCorrelations(AcqFiles=[], ax=None):
    colors = ['b', 'g', 'r', 'y']
    if ax != None:
        ax[0].clear()
    linex = numpy.array([numpy.ones(50), numpy.linspace(0, 0.3)])
    difference = {}
    strehl = {}
    for i in range(4):
        A = []
        x = []
        strehl[i] = []
        for Acq in AcqFiles:
            for SR, flux, dx, dy, sc in zip(Acq.newStrehl.data[i], Acq.TOTALFLUX_FT.data[i], Acq.newSC_FIBER_DX.data[i], Acq.newSC_FIBER_DY.data[i], Acq.newSCALE.data[i]):
                if not(numpy.isnan(SR)):
                    strehl[i].append(SR)
                    #A.append([1.0, SR/numpy.max(Acq.newStrehl.data[i])])
                    #x.append(flux / Acq.TOTALFLUX_FT.median[i])
                    A.append([1.0, SR])
                    #A.append([1.0, SR, (dx**2.0+dy**2.0)**0.5])
                    #A.append([1.0, SR, dx, dy])
                    #A.append([1.0, SR, dx, dy, sc])
                    x.append(flux )
                    if ax != None:
                        #ax.scatter([SR/numpy.max(Acq.newStrehl.data[i])], [flux/Acq.TOTALFLUX_FT.median[i]], c=colors[i])
                        ax[0].scatter([SR], [flux], c=colors[i])

        A = numpy.array(A)
        x = numpy.array(x)
        B = numpy.linalg.pinv(A)
        fit = B.dot(x)
        if ax != None:
            #ax.plot(linex[1], fit.dot(linex), c=colors[i])
            ax[0].scatter(A[:,1], fit.dot(A.T), c=colors[i], marker='.')

            for blah, junk, fart in zip(A, x, fit.dot(A.T)):
                ax[0].plot([blah[1], blah[1]], [junk, fart], c=colors[i])
        difference[i] = numpy.mean(numpy.sqrt(((x - fit.dot(A.T))**2.0/fit.dot(A.T))))
        print difference[i]
        strehl[i] = numpy.array(strehl[i])

    for i in range(4):
        for j in range(i+1, 4):
            print "Correlation between Strehls of Telescopes ", i, " and ", j, ": ", numpy.corrcoef(strehl[i], strehl[j])

    if ax != None:
        ax[1].clear()
        ax[1].plot(strehl[0], c=colors[0])
        ax[1].plot(strehl[1], c=colors[1])
        ax[1].plot(strehl[2], c=colors[2])
        ax[1].plot(strehl[3], c=colors[3])

    


fig1 = pyplot.figure(0)
fig1.clear()
ax1 = fig1.add_axes([0.1, 0.1, 0.8, 0.8])
ax1.set_title("Science Camera")
fig2 = pyplot.figure(1)
fig2.clear()
ax2 = fig2.add_axes([0.1, 0.1, 0.8, 0.8])
ax2.set_title("BCI Flux")
fig3 = pyplot.figure(2)
fig3.clear()
ax3 = fig3.add_axes([0.1, 0.1, 0.8, 0.8])
ax3.set_title("Plate Scale")
fig4 = pyplot.figure(3)
fig4.clear()
ax4 = fig4.add_axes([0.1, 0.1, 0.8, 0.8])
ax4.set_title("Fiber")
fig5 = pyplot.figure(4)
fig5.clear()
ax5 = fig5.add_axes([0.1, 0.1, 0.8, 0.8])
ax5.set_title("Strehl")

axes = [ax1, ax2, ax3, ax4, ax5]

files = glob.glob('/home/cdeen/Data/GRAVITY/2017-08/2017-08-09/reduced/GRAVI*dualscip2vmred.fits')

AcqFiles = []

for df in files:
    fiberData = pyfits.getdata(df, ext=16)
    fluxData = pyfits.getdata(df, ext=8)
    AcqFiles.append(AcqCamData(fiberData=fiberData, fluxData=fluxData))
    AcqFiles[-1].binData()
    blah = AcqFiles[-1].calcMedian()


findCorrelations(AcqFiles, ax=[ax1, ax2])

#Acq.binData()
#Acq.plot(axes=axes)
#Acq.findCorrelations(ax1)

fig1.show()
#fig2.show()
#fig3.show()
#fig4.show()
#fig5.show()
