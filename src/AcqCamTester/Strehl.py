import astropy.io.fits as pyfits
from matplotlib import pyplot
import numpy
import scipy
import glob
import sys

sys.path.append('../')
import CIAO_DatabaseTools
import Graffity
import astropy.time as aptime

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
    def __init__(self, fiberData=None, fluxData=None, FTMag=0.0, AcqDit=0.0,
            CIAO_Data = None):
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

        self.FTMag = FTMag
        self.AcqDit = AcqDit
        self.CIAO_Data = CIAO_Data

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
                    ax.plotlways (self.BCI_Time.data[i], self.TOTALFLUX_SC.data[i]/medianTOTALFLUX_SC[i], c=color)
                    ax.plot(self.BCI_Time.data[i], self.TOTALFLUX_FT.data[i]/medianTOTALFLUX_FT[i], c=color)
                elif name == 'SCALE':
                    ax.errorbar(self.BCI_Time.data[i], self.newSCALE.data[i]-medianSCALE[i], yerr=self.newSCALE.error[i], c=color)
                elif name == 'FIBER':
                    ax.errorbar(self.BCI_Time.data[i], self.newSC_FIBER_DX.data[i], yerr=self.newSC_FIBER_DX.error[i], c=color, ls='-')
                    ax.errorbar(self.BCI_Time.data[i], self.newSC_FIBER_DY.data[i], yerr=self.newSC_FIBER_DY.error[i], c=color, ls='-.')
                elif name == 'STREHL':
                    ax.plot(self.BCI_Time.data[i], self.newStrehl.data[i]/medianStrehl[i], c=color)


def findCorrelations(AcqFiles=[], ax=None):
    colors = ['b', 'g', 'r', 'y']
    if ax != None:
        ax[0].clear()
    linex = numpy.array([numpy.ones(50), numpy.linspace(0, 0.3)])
    difference = {}
    strehl = {}
    cov = {}

    for i,j in zip(range(4), [3, 2, 1,0]):
        A = []
        x = []
        strehl[i] = []
        for Acq in AcqFiles:
            for SR, flux, dx, dy, sc in zip(Acq.newStrehl.data[i],
                    Acq.TOTALFLUX_FT.data[i], Acq.newSC_FIBER_DX.data[i],
                    Acq.newSC_FIBER_DY.data[i], Acq.newSCALE.data[i]):
                strehlAccumulator = []
                if (not(numpy.isnan(SR)) and (SR > 0) and (SR < 0.5) and
                        (Acq.FTMag == 9.7)):
                    if (flux > -150000):
                        strehl[i].append(SR)
                        strehlAccumulator.append(SR)
                        #A.append([1.0, SR/numpy.max(Acq.newStrehl.data[i])])
                        #x.append(flux / Acq.TOTALFLUX_FT.median[i])
                        #A.append([1.0, SR])
                        #A.append([1.0, SR, 10.0**(-2.5*Acq.FTMag)*Acq.AcqDit])
                        A.append([1.0, SR, 10.0**(-2.5*Acq.FTMag)*Acq.AcqDit,
                            Acq.CIAO_Data[j]["Strehl"],
                            Acq.CIAO_Data[j]["Seeing"],
                            Acq.CIAO_Data[j]["Tau0"]])
                            #Acq.CIAO_Data[j]["WindSp"]])
                        #A.append([1.0, SR, (dx**2.0+dy**2.0)**0.5])
                        #A.append([1.0, SR, dx, dy])
                        #A.append([1.0, SR, dx, dy, sc])
                        x.append(flux )
                        if ax != None:
                            #ax.scatter([SR/numpy.max(Acq.newStrehl.data[i])], [flux/Acq.TOTALFLUX_FT.median[i]], c=colors[i])
                            ax[0].scatter([SR], [flux], c=colors[i])
            #if ax!= None:
            #    ax[0].scatter([numpy.mean(numpy.array(strehlAccumulator))],
            #            [Acq.CIAO_Data[j]["Strehl"]], c=colors[i])
                            

        """
        #A = numpy.array(A, dtype=float)
        #x = numpy.array(x)
        #B = numpy.linalg.pinv(A)
        #fit = B.dot(x)
        #if ax != None:
            #ax.plot(linex[1], fit.dot(linex), c=colors[i])
            #ax[0].scatter(A[:,1], fit.dot(A.T), c=colors[i], marker='.')

            #for blah, junk, fart in zip(A, x, fit.dot(A.T)):
            #    ax[0].plot([blah[1], blah[1]], [junk, fart], c=colors[i])
        #difference[i] = numpy.mean(numpy.sqrt(((x - fit.dot(A.T))**2.0/fit.dot(A.T))))
        #print difference[i]
        #print fit, A.shape

        #n = A.shape[0]
        #avg = numpy.mean(A, axis=0)
        #cov[i] = 1.0/n * (A - avg).T.dot(A-avg)
        #strehl[i] = numpy.array(strehl[i])

    ax[0].set_xlabel("H-band AcqCam Strehl")
    ax[0].set_ylabel("FT Flux")
    #for i in range(4):
    #    for j in range(i+1, 4):
    #        print "Correlation between Strehls of Telescopes ", i, " and ", j, ": ", numpy.corrcoef(strehl[i], strehl[j])
    """
    if ax != None:
        ax[1].clear()
        ax[1].plot(strehl[0], c=colors[0])
        ax[1].plot(strehl[1], c=colors[1])
        ax[1].plot(strehl[2], c=colors[2])
        ax[1].plot(strehl[3], c=colors[3])

    return 0#cov
    #"""


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

fig6 = pyplot.figure(5)
fig6.clear()
ax61 = fig6.add_axes([0.1, 0.1, 0.4, 0.4])
ax62 = fig6.add_axes([0.1, 0.5, 0.4, 0.4])
ax63 = fig6.add_axes([0.5, 0.1, 0.4, 0.4])
ax64 = fig6.add_axes([0.5, 0.5, 0.4, 0.4])
fig6.suptitle("Covariances")

axes = [ax1, ax2, ax3, ax4, ax5]

#files = glob.glob('/home/grav/cdeen/GRAVITY/reduced/GRAVI*singlescip2vmred.fits')
files = glob.glob('/home/grav/cdeen/GRAVITY/reduced/GRAVI.2017-07*dualscip2vmred.fits')

AcqFiles = []

CDB = CIAO_DatabaseTools.CIAO_Database()
CIAO_DataLoggers = CDB.query(keywords = ["TAU0", "STREHL", "SEEING", "WINDSP"])

#for df in['/home/grav/cdeen/GRAVITY/reduced/GRAVI.2017-08-09T01:12:28.341_dualscip2vmred.fits']:
for df in files:
    data = pyfits.open(df)
    if len(data) > 16:
        fiberData = pyfits.getdata(df, ext=16)
        fluxData = pyfits.getdata(df, ext=10)
        FTMag = pyfits.getheader(df).get('ESO FT ROBJ MAG')
        AcqDit = data[0].header.get('ESO DET1 SEQ1 DIT')
        print FTMag, AcqDit
        time = aptime.Time(data[0].header.get('ESO PCR ACQ START'))
        CIAO_Data = {}
        good = True
        for i in [1,2, 3, 4]:
            CIAO_Data[i-1] = {}
            times = numpy.array(CIAO_DataLoggers[i][:,-4], dtype=float)
            closest = numpy.argsort(numpy.abs(times-time.mjd))[0]
            if numpy.abs(times[closest] - time.mjd) > (1.0/(24.0*3600/30.0)):
                print "Ruh-roh!  Difference in time is greater than 30 seconds!"
                good = False
            else:
                print "Time Difference = ", (times[closest] - time.mjd)*24*3600.0
                CIAO_Data[i-1]["Tau0"] = CIAO_DataLoggers[i][closest,0]
                CIAO_Data[i-1]["Strehl"] = CIAO_DataLoggers[i][closest,1]
                CIAO_Data[i-1]["Seeing"] = CIAO_DataLoggers[i][closest,2]
                CIAO_Data[i-1]["WindSp"] = CIAO_DataLoggers[i][closest,3]
                print CIAO_Data[i-1]
        if good:
            AcqFiles.append(AcqCamData(fiberData=fiberData,
                 fluxData=fluxData,FTMag=FTMag, AcqDit=AcqDit, CIAO_Data = CIAO_Data))
            AcqFiles[-1].binData()
            blah = AcqFiles[-1].calcMedian()
    else:
        print "%s does not have all required tables!" % df
    data.close()


covar = findCorrelations(AcqFiles, ax=[ax1, ax2])

#for cov, ax in zip(covar.keys(), [ax61, ax62, ax63, ax64]):
#    ax.matshow(numpy.log10(covar[cov][1:,1:]**2.0), vmin=-20)

#Acq.binData()
#Acq.plot(axes=axes)
#Acq.findCorrelations(ax1)

fig1.show()
fig2.show()
#fig3.show()
#fig4.show()
#fig5.show()
#fig6.show()
