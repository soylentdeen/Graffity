import numpy
import scipy
import astropy.io.fits as pyfits
import os
from scipy.misc import factorial as fac
import scipy.interpolate as interp
import scipy.fftpack as fftpack
import matplotlib.pyplot as pyplot
from scipy import optimize
from scipy.ndimage import rotate
from PIL import Image
from scipy import signal
from scipy import linalg
from scipy.special import gamma
from scipy.ndimage.measurements import center_of_mass
from scipy.ndimage.filters import maximum_filter
from scipy.ndimage.morphology import generate_binary_structure, binary_erosion
import time
from datetime import datetime
from astropy import time as aptime

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
        else:
            for i in range(4):
                self.data[i] = numpy.array([])

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

    def append(self, other):
        for i in self.data.keys():
            self.data[i] = numpy.append(self.data[i], other.data[i])
            if self.error != None:
                self.error[i] = numpy.append(self.error[i], other.error[i])
        return self

    def median(self):
        for i in self.data.keys():
            self.data[i] = numpy.median(self.data[i])

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
            #time.data[i] = time.data[i][time.data[i] != 0]
            if self.error != None:
                new.error[i] = []
            startTime = otherTime.data[i][0] - numpy.mean(numpy.diff(otherTime.data[i]))
            stopTime = (otherTime.data[i][1] + otherTime.data[i][0]) / 2.0
            included = (startTime < time.data[i]) & (stopTime > time.data[i]) & (self.data[i] != 0)
            if self.error != None:
                new.data[i].append(numpy.sum(self.data[i][included]/self.error[i][included])/numpy.sum(
                              1.0/self.error[i][included]))
                new.error[i].append(numpy.sum(self.error[i][included]**-2.0)**-0.5)
            else:
                new.data[i].append(numpy.mean(self.data[i][included]))
            for t in range(len(otherTime.data[i])-2):
                startTime = (otherTime.data[i][t]+otherTime.data[i][t+1])/2.0
                stopTime = (otherTime.data[i][t+1]+otherTime.data[i][t+2])/2.0
                included = (startTime < time.data[i]) & (stopTime > time.data[i]) & (self.data[i] != 0)
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
            included = (startTime < time.data[i]) & (stopTime > time.data[i]) & (self.data[i] != 0)
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
    def __init__(self, fiberData=None, fluxData=None, FTMag=0.0, SCMag= 0.0, AcqDit=0.0,
            CIAO_Data =None, imagingData=None):
        self.Time = BCI_Data(data=fiberData.field('TIME'))
        self.Strehl = BCI_Data(data=fiberData.field('FIELD_STREHL'))
        self.SC_X = BCI_Data(data=fiberData.field('FIELD_SC_X'), error=fiberData.field('FIELD_SC_XERR'))
        self.SC_Y = BCI_Data(data=fiberData.field('FIELD_SC_Y'), error=fiberData.field('FIELD_SC_YERR'))
        self.FT_X = BCI_Data(data=fiberData.field('FIELD_FT_X'), error=fiberData.field('FIELD_FT_XERR'))
        self.FT_Y = BCI_Data(data=fiberData.field('FIELD_FT_Y'), error=fiberData.field('FIELD_FT_YERR'))
        self.SCALE = BCI_Data(data=fiberData.field('FIELD_SCALE'), error=fiberData.field('FIELD_SCALEERR'))
        self.SC_FIBER_DX = BCI_Data(data=fiberData.field('FIELD_FIBER_DX'), error=fiberData.field('FIELD_FIBER_DXERR'))
        self.SC_FIBER_DY = BCI_Data(data=fiberData.field('FIELD_FIBER_DY'), error=fiberData.field('FIELD_FIBER_DYERR'))
        self.Pupil_X = BCI_Data(data=fiberData.field('PUPIL_X'))
        self.Pupil_Y = BCI_Data(data=fiberData.field('PUPIL_Y'))
        self.Pupil_Z = BCI_Data(data=fiberData.field('PUPIL_Z'))
        self.Pupil_R = BCI_Data(data=fiberData.field('PUPIL_R'))
        self.Pupil_U = BCI_Data(data=fiberData.field('PUPIL_U'))
        self.Pupil_V = BCI_Data(data=fiberData.field('PUPIL_V'))
        self.Pupil_W = BCI_Data(data=fiberData.field('PUPIL_W'))
        self.OPD_Pupil = BCI_Data(data=fiberData.field('OPD_PUPIL'))

        self.BCI_Time = BCI_Data(fluxData['TIME'])
        self.TOTALFLUX_SC = BCI_Data(fluxData['TOTALFLUX_SC'])
        self.TOTALFLUX_FT = BCI_Data(fluxData['TOTALFLUX_FT'])
        self.FDDL = BCI_Data(fluxData['FDDL'])
        self.FT_Pos = BCI_Data(fluxData['FT_POS'])
        self.SC_Pos = BCI_Data(fluxData['SC_POS'])

        self.FTMag = FTMag
        self.SCMag = SCMag
        self.AcqDit = AcqDit
        self.CIAO_Data = CIAO_Data
        self.AcqFT = None
        self.AcqSC = None
        if imagingData != None:
            #self.xi, self.yi = numpy.meshgrid(numpy.arange(20),
            #        numpy.arange(20))
            #self.calcEllipse(imagingData)
            self.calcMaxPixel(imagingData)

    def binData(self):
        self.newStrehl = self.Strehl.rebin(self.Time, self.BCI_Time)
        self.newSC_X = self.SC_X.rebin(self.Time, self.BCI_Time)
        self.newSC_Y = self.SC_Y.rebin(self.Time, self.BCI_Time)
        self.newFT_X = self.FT_X.rebin(self.Time, self.BCI_Time)
        self.newFT_Y = self.FT_Y.rebin(self.Time, self.BCI_Time)
        self.newSCALE = self.SCALE.rebin(self.Time, self.BCI_Time)
        self.newSC_FIBER_DX = self.SC_FIBER_DX.rebin(self.Time, self.BCI_Time)
        self.newSC_FIBER_DY = self.SC_FIBER_DY.rebin(self.Time, self.BCI_Time)
        self.newPupilX = self.Pupil_X.rebin(self.Time, self.BCI_Time)
        self.newPupilY = self.Pupil_Y.rebin(self.Time, self.BCI_Time)
        self.newPupilZ = self.Pupil_Z.rebin(self.Time, self.BCI_Time)
        self.newPupilR = self.Pupil_R.rebin(self.Time, self.BCI_Time)
        self.newPupilU = self.Pupil_U.rebin(self.Time, self.BCI_Time)
        self.newPupilV = self.Pupil_V.rebin(self.Time, self.BCI_Time)
        self.newPupilW = self.Pupil_W.rebin(self.Time, self.BCI_Time)
        self.newOPDPupil = self.OPD_Pupil.rebin(self.Time, self.BCI_Time)
        if self.AcqFT != None:
            self.newAcqFT = self.AcqFT.rebin(self.Time, self.BCI_Time)
        if self.AcqSC != None:
            self.newAcqSC = self.AcqSC.rebin(self.Time, self.BCI_Time)


    def calcMedian(self):
        medianTOTALFLUX_SC = self.TOTALFLUX_SC.getMedian()
        medianTOTALFLUX_FT = self.TOTALFLUX_FT.getMedian()
        medianSCALE = self.SCALE.getMedian()
        medianSC_FIBER_DX = self.SC_FIBER_DX.getMedian()
        medianSC_FIBER_DY = self.SC_FIBER_DY.getMedian()
        medianStrehl = self.Strehl.getMedian()

    def newImg(self, p):
        return p[0]*numpy.exp(-(self.xi-p[1])**2.0/p[2] -
                ((self.yi-p[3])**2.0/p[4]))

    def makeEllipse(self, p, img):
        return (self.newImg(p).flatten() - img.flatten())**2.0

    def fitEllipse(self, img):
        img /= numpy.max(img)
        c = 1.0
        rx = 1.0
        ry = 1.0
        xc = 10.0
        yc = 10.0
        fit, result = optimize.leastsq(self.makeEllipse, [c, xc, rx, yc, ry], args=(img))
        ratio = fit[2]/fit[4]
        return ratio, fit[1], fit[3]

    def calcMaxPixel(self, imagingData):
        header = imagingData[0]
        images = imagingData[1]
        subImageSize = 250
        FT = []
        SC = []
        for j, im in enumerate(images):
            for i in [1, 2, 3, 4]:
                startX = header.get('ESO DET1 FRAM%d STRX' %i)
                startY = header.get('ESO DET1 FRAM%d STRY' %i)
                xcoord = (self.FT_X.data[i-1][j] - startX) + (i-1)*subImageSize
                ycoord = self.FT_Y.data[i-1][j] - startY
                try:
                    FT.append(numpy.max(im[int(ycoord-3):int(ycoord+3),
                        int(xcoord-3):int(xcoord+3)]))
                except:
                    FT.append(numpy.nan)
                xcoord = (self.SC_X.data[i-1][j] - startX) + (i-1)*subImageSize
                ycoord = self.SC_Y.data[i-1][j] - startY
                try:
                    SC.append(numpy.max(im[int(ycoord-3):int(ycoord+3),
                         int(xcoord-3):int(xcoord+3)]))
                except:
                    SC.append(numpy.nan)
        self.AcqFT = BCI_Data(data = numpy.array(FT))
        self.AcqSC = BCI_Data(data = numpy.array(SC))
        
    def calcEllipse(self, imagingData):
        header = imagingData[0]
        images = imagingData[1]
        subImageSize = 250
        ellipse = []
        xpos = []
        ypos = []
        for j, im in enumerate(images):
            for i in [1, 2, 3, 4]:
                startX = header.get('ESO DET1 FRAM%d STRX' %i)
                startY = header.get('ESO DET1 FRAM%d STRY' %i)
                xcoord = (self.FT_X.data[i-1][j] - startX) + (i-1)*subImageSize
                ycoord = self.FT_Y.data[i-1][j] - startY
                FT_Fit = self.fitEllipse(im[ycoord-10:ycoord+10,
                    xcoord-10:xcoord+10])
                ellipse.append(FT_Fit[0])
                xpos.append(FT_Fit[1]-10 + xcoord)
                ypos.append(FT_Fit[2]-10 + ycoord)

            del(im)

        self.Ellipse=BCI_Data(data=ellipse)
        self.Fit_X = BCI_Data(data=xpos)
        self.Fit_Y = BCI_Data(data=ypos)
        

class PSF( object ):
    def __init__(self, sizeInPix=20, lam=1.6, dlam=0.3, pscale=17.0, M1=8.0, M2=1.3, nLambdaSteps=9):
        self.sizeInPix = sizeInPix
        self.lam = lam*1.0e-6
        self.dlam = dlam*1.0e-6
        self.pscale = pscale * 2.0 *numpy.pi/360.0/60.0**2.0/1.0e3
        #self.pscale = pscale/360.0/60.0**2.0/1.0e3
        self.M1 = M1
        self.M2 = M2
        self.nLambdaSteps = nLambdaSteps
        self.fmax = self.M1 * self.pscale * self.sizeInPix/self.lam
        self.Zernikes = []
        for i in range(50):
            self.Zernikes.append(zernikeMode(i, 0.0))

    def setZern(self, zern):
        """
            Sets the magnitudes of the Zernike components.
        """
        nzern = len(zern)
        for i in range(nzern):
            self.Zernikes[i].setMag(zern[i]*2.0*numpy.pi/self.lam)

    def Sinc(self, x):
        if numpy.abs(x) < 1e-4:
            return 1.0
        else:
            return numpy.sin(x)/x

    def H1(self, f, u, v):
        if numpy.abs(1.0-v) < 1e-12:
            e = 1
        else:
            e = -1
        return v**2.0/numpy.pi *numpy.arccos((f/v) * (1.0+e*(1.0-u**2.0)/(4*f**2.0)))

    def H2(self, f, u):
        a = 2.0*f/(1.0+u)
        b = (1.0-u)/(2.0*f)
        return -1.0*(f/numpy.pi)*(1.0+u)*numpy.sqrt((1.0-a**2.0)*(1.0-b**2.0))

    def G(self, f, u):
        if f <= (1.0 - u)/2.0:
            return u**2.0
        elif f >= (1.0 + u)/2.0:
            return 0.0
        else:
            return self.H1(f, u, 1.0) + self.H1(f, u, u) + self.H2(f, u)

    def TelOTF(self, f, u):
        retval = (self.G(f, 1.0) + u**2.0*self.G(f/u, 1) - 2*self.G(f, u))/(1.0 - u**2.0)
        return retval

    def generateOTF(self, ax=None):
        total = numpy.zeros([self.sizeInPix, self.sizeInPix], dtype=numpy.complex)

        for i in range(self.sizeInPix):
            for j in range(self.sizeInPix):
                y = j - self.sizeInPix/2.0 + 0.5
                x = i - self.sizeInPix/2.0 + 0.5
                r = numpy.sqrt(x**2.0 + y**2.0)
                value = 0.0
                phase = 0.0
                for k in numpy.arange(self.nLambdaSteps):
                    l = self.lam - self.dlam*(k - self.nLambdaSteps/2.0)/(self.nLambdaSteps-1)
                    fc = self.fmax*self.lam/l
                    f = r/fc
                    if f < 1:
                        if r < 0.1:
                            value += 1.0/self.nLambdaSteps
                        else:
                            value += (self.TelOTF(f, self.M2/self.M1)*
                                            self.Sinc(numpy.pi*x/self.sizeInPix)*
                                            self.Sinc(numpy.pi*y/self.sizeInPix))/self.nLambdaSteps
                        if k==0:
                            phi = numpy.arctan2(y, x)
                            for z in self.Zernikes:
                                phase += z.zernike(f, phi)

                #total[i, j] = value*numpy.exp(numpy.complex(0.0, 2.0*numpy.pi*phase))
                total[i, j] = value*numpy.exp(numpy.complex(0.0, phase))
        self.OTF = total
        self.PSF = numpy.fft.fftshift(numpy.fft.fft2(self.OTF))
        self.PSF = ((self.PSF*self.PSF.conj())**0.5).real
        self.PSF = self.normalize(self.PSF)
        if ax != None:
            ax.matshow(numpy.imag(self.OTF))
            ax.figure.show()
            raw_input()

    def getPSF(self):
        return self.PSF

    def normalize(self, image):
        return image/numpy.sum(image)

    def calcStrehl(self, cutout):
        factor = numpy.max(cutout) / numpy.max(self.PSF)
        if factor > 1.0:
            raise
        return factor
        


class SubAperture( object ):
    def __init__(self, index):
        self.index = index
        self.x = []
        self.y = []
        self.z = []

    def addFluxPoint(self, tip, tilt, flux):
        self.x.append(tip)
        self.y.append(tilt)
        self.z.append(flux)

    def interpolate(self, tipRange=[-1, 1], tiltRange=[-1,1], spacing=0.01):
        interpol = interp.interp2d(self.x, self.y, self.z, fill_value=0.0)
        gridx = numpy.arange(tipRange[0], tipRange[1], spacing)
        gridy = numpy.arange(tiltRange[0], tiltRange[1], spacing)

        return interpol(gridx, gridy)


class FOVFrame( object ):
    def __init__(self, naps = 68):
        self.subapertures = []
        for i in range(naps):
            self.subapertures.append(SubAperture(i))

    def addFrame(self, tip, tilt, frame):
        for subap, flux in zip(self.subapertures, frame):
            subap.addFluxPoint(tip, tilt, flux)

    def publish(self, ax):
        for subap in self.subapertures:
            postageStamp = subap.interpolate(tipRange=[-1, 1], 
                                tiltRange=[-1, 1], spacing=0.01)
            ax.clear()
            ax.imshow(postageStamp)
            ax.figure.show()
            raw_input()
        

class LaserBeaconData( object ):
    def __init__(self, df=''):
        self.df = df
        self.header = pyfits.getheader(df)
        self.HDU = pyfits.open(df)
        self.Time = self.HDU[5].data['TIME']/1e6
        self.Piezo1 = self.HDU[5].data['Piezo1']
        self.Piezo2 = self.HDU[5].data['Piezo2']
        self.Piezo3 = self.HDU[5].data['Piezo3']
        self.Piezo4 = self.HDU[5].data['Piezo4']
        self.PSD1 = self.HDU[5].data['PSD1']
        self.PSD2 = self.HDU[5].data['PSD2']
        self.PSD3 = self.HDU[5].data['PSD3']
        self.PSD4 = self.HDU[5].data['PSD4']

        self.loopRate = 1.0/numpy.mean(numpy.diff(self.Time))
        self.R = 0.5

    def computePSDs(self):
        NWindows = computeNWindows(self.Piezo1, self.loopRate, self.R)
        self.Piezo1 = FourierDetrend(self.Piezo1, 'TwoPoints')

        WindowSize = 20.0
        N = int(WindowSize*self.loopRate)

        self.PiezoPower = {1:[], 2:[], 3:[], 4:[]}
        self.PSDPower = {1:[], 2:[], 3:[], 4:[]}
        i = 0
        while i+N < len(self.Time):
            time = self.Time[i:i+N]
            
            PiezoPower = {}
            PSDPower = {}
            Piezo1 = self.Piezo1[i:i+N,:]
            Piezo2 = self.Piezo2[i:i+N,:]
            Piezo3 = self.Piezo3[i:i+N,:]
            Piezo4 = self.Piezo4[i:i+N,:]
            PSD1 = self.PSD1[i:i+N,:]
            PSD2 = self.PSD2[i:i+N,:]
            PSD3 = self.PSD3[i:i+N,:]
            PSD4 = self.PSD4[i:i+N,:]
            i = i+N
            NFrames = Piezo1.shape[0]
            PiezoPower[1] = numpy.abs(numpy.fft.fft(Piezo1.T/NFrames))**2.0
            PiezoPower[2] = numpy.abs(numpy.fft.fft(Piezo2.T/NFrames))**2.0
            PiezoPower[3] = numpy.abs(numpy.fft.fft(Piezo3.T/NFrames))**2.0
            PiezoPower[4] = numpy.abs(numpy.fft.fft(Piezo4.T/NFrames))**2.0
            PSDPower[1] = numpy.abs(numpy.fft.fft(PSD1.T/NFrames))**2.0
            PSDPower[2] = numpy.abs(numpy.fft.fft(PSD2.T/NFrames))**2.0
            PSDPower[3] = numpy.abs(numpy.fft.fft(PSD3.T/NFrames))**2.0
            PSDPower[4] = numpy.abs(numpy.fft.fft(PSD4.T/NFrames))**2.0
            Freq = numpy.fft.fftfreq(NFrames, d=1.0/self.loopRate)

            PiezoPower[1] = numpy.fft.fftshift(PiezoPower[1], axes=1)
            PiezoPower[2] = numpy.fft.fftshift(PiezoPower[2], axes=1)
            PiezoPower[3] = numpy.fft.fftshift(PiezoPower[3], axes=1)
            PiezoPower[4] = numpy.fft.fftshift(PiezoPower[4], axes=1)
            PSDPower[1] = numpy.fft.fftshift(PSDPower[1], axes=1)
            PSDPower[2] = numpy.fft.fftshift(PSDPower[2], axes=1)
            PSDPower[3] = numpy.fft.fftshift(PSDPower[3], axes=1)
            PSDPower[4] = numpy.fft.fftshift(PSDPower[4], axes=1)
            Freq = numpy.fft.fftshift(Freq)
            PiezoPower[1] = 2.0*PiezoPower[1][:,Freq >= 0]
            PiezoPower[2] = 2.0*PiezoPower[2][:,Freq >= 0]
            PiezoPower[3] = 2.0*PiezoPower[3][:,Freq >= 0]
            PiezoPower[4] = 2.0*PiezoPower[4][:,Freq >= 0]
            PSDPower[1] = 2.0*PSDPower[1][:,Freq >= 0]
            PSDPower[2] = 2.0*PSDPower[2][:,Freq >= 0]
            PSDPower[3] = 2.0*PSDPower[3][:,Freq >= 0]
            PSDPower[4] = 2.0*PSDPower[4][:,Freq >= 0]
            Freq = Freq[Freq >= 0]
            PiezoPower[1][:,0] /= 2.0
            PiezoPower[2][:,0] /= 2.0
            PiezoPower[3][:,0] /= 2.0
            PiezoPower[4][:,0] /= 2.0
            PSDPower[1][:,0] /= 2.0
            PSDPower[2][:,0] /= 2.0
            PSDPower[3][:,0] /= 2.0
            PSDPower[4][:,0] /= 2.0
            PiezoPower[1] = PiezoPower[1] / numpy.mean(numpy.diff(Freq))
            PiezoPower[2] = PiezoPower[2] / numpy.mean(numpy.diff(Freq))
            PiezoPower[3] = PiezoPower[3] / numpy.mean(numpy.diff(Freq))
            PiezoPower[4] = PiezoPower[4] / numpy.mean(numpy.diff(Freq))
            PSDPower[1] = PSDPower[1] / numpy.mean(numpy.diff(Freq))
            PSDPower[2] = PSDPower[2] / numpy.mean(numpy.diff(Freq))
            PSDPower[3] = PSDPower[3] / numpy.mean(numpy.diff(Freq))
            PSDPower[4] = PSDPower[4] / numpy.mean(numpy.diff(Freq))

            NSamplesPerWindow = numpy.floor(PiezoPower[1].shape[1]/NWindows)
            TotalSampleNumber = NWindows * NSamplesPerWindow
            self.PiezoPower[1].append(PiezoPower[1].T)
            self.PiezoPower[2].append(PiezoPower[2].T)
            self.PiezoPower[3].append(PiezoPower[3].T)
            self.PiezoPower[4].append(PiezoPower[4].T)
            self.PSDPower[4].append(PSDPower[4].T)
            self.PSDPower[3].append(PSDPower[3].T)
            self.PSDPower[2].append(PSDPower[2].T)
            self.PSDPower[1].append(PSDPower[1].T)

        self.Freq = Freq
        for i in [1, 2, 3, 4]:
            self.PiezoPower[i] = numpy.mean(numpy.array(self.PiezoPower[i]), axis=0)
            self.PSDPower[i] = numpy.mean(numpy.array(self.PSDPower[i]), axis=0)
        



class PixelBuffer( object ):
    def __init__(self, df=''):
        self.df = df
        self.directory = os.path.dirname(self.df)
        self.header = pyfits.getheader(df)
        self.data = pyfits.getdata(columns)
        self.FrameCounter = self.data.field('FrameCounter')
        self.Pixels = self.data.field('Pixels')

    def computeCentroids(self):
        centroids = []
        for pix in self.Pixels:
            centroids.append(self.computeCentroids(pix))

        self.centroids = numpy.array(centroids)

class WFS_Frame( object ):
    def __init__(self, intensities = None, gradients = None):
        self.intensities = intensities
        self.gradients = gradients
        self.subap = numpy.array([
                [False, False, True, True, True, True, True, False, False],
                [False, True, True, True, True, True, True, True, False],
                [True, True, True, True, True, True, True, True, True],
                [True, True, True, True, True, True, True, True, True],
                [True, True, True, True, False, True, True, True, True],
                [True, True, True, True, True, True, True, True, True],
                [True, True, True, True, True, True, True, True, True],
                [False, True, True, True, True, True, True, True, False],
                [False, False, True, True, True, True, True, False, False]])
        self.image = numpy.zeros(self.subap.shape)
        self.generateMasks()
        
    def generateMasks(self):
        self.innerRing = numpy.ones(self.subap.shape) == 0.0
        self.outerRing = numpy.ones(self.subap.shape) == 0.0

        self.innerRing[3:6,3:6] = True
        self.outerRing[0, 2:7] = True
        self.outerRing[-1, 2:7] = True
        self.outerRing[2:7, 0] = True
        self.outerRing[2:7, -1] = True

    def generateIntensityImage(self, intensities = None):
        if intensities != None:
            self.intensities = intensities
        self.image[self.subap] = self.intensities
        self.innerRingImage = numpy.zeros(self.subap.shape)
        self.innerRingImage[self.innerRing] = self.image[self.innerRing]
        self.outerRingImage = numpy.zeros(self.subap.shape)
        self.outerRingImage[self.outerRing] = self.image[self.outerRing]


    def plotIntensities(self, ax=None):
        if ax == None:
            fig = pyplot.figure(0)
            ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        ax.matshow(self.image)
        ax.figure.show()

    def plotGradients(self, ax=None):
        if ax == None:
            fig = pyplot.figure(0)
            fig.clear()
            ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
            ax.clear()
        subapwidth = 8
        counter = 0
        for row in range(len(self.subap)):
            for col in range(len(self.subap[0])):
                if self.subap[row][col]:
                    centerx = row*subapwidth+subapwidth/2.0-0.5
                    centery = col*subapwidth+subapwidth/2.0-0.5
                    distx = centerx - self.gradients[counter][1]
                    disty = centery - self.gradients[counter][0]
                    ax.plot([centery, centery+disty*100],[centerx,
                             centerx+distx*100])
                    counter += 1
        print self.gradients[5][1]
        fig.show()

class GRAVITY_Dual_P2VM( object ):
    def __init__(self, fileBase = '', startTime=0.0, CIAO_Data=None,
            processAcqCamData=False):
        return self

    def getPupilMotion(self):
        self.TIME = (self.data.field('TIME')[::4]*1e-6)/86400.0+self.startTime.mjd
        self.PUPIL_NSPOT = self.data.field('PUPIL_NSPOT')[::4]
        self.PUPIL_X = {}
        self.PUPIL_Y = {}
        self.PUPIL_Z = {}
        self.PUPIL_R = {}
        self.PUPIL_U = {}
        self.PUPIL_V = {}
        self.PUPIL_W = {}
        
        for UT in [1, 2, 3, 4]:
            self.PUPIL_X[UT] = self.data.field('PUPIL_X')[UT-1::4]
            self.PUPIL_Y[UT] = self.data.field('PUPIL_Y')[UT-1::4]
            self.PUPIL_Z[UT] = self.data.field('PUPIL_Z')[UT-1::4]
            self.PUPIL_R[UT] = self.data.field('PUPIL_R')[UT-1::4]
            self.PUPIL_U[UT] = self.data.field('PUPIL_U')[UT-1::4]
            self.PUPIL_V[UT] = self.data.field('PUPIL_V')[UT-1::4]
            self.PUPIL_W[UT] = self.data.field('PUPIL_W')[UT-1::4]

    def computeOPDPeriodograms(self):
        M_matrix = numpy.array([-1.,1.,0.0,0.0,-1.,0.0,1.,0.0,-1.,0.0,0.0,1.,0.0,-1.,1.,0.0,0.0,-1.,0.0,1.,0.0,0.0,-1.,1.]);
        M_matrix = M_matrix.reshape((6, 4))

        opdc_kalman_pizeo_opd = numpy.dot(M_matrix,
                self.opdc_kalman_piezo.T).T
        PSD_k = [1, 2, 3, 4, 5, 6]
        for baseline in range(0, 6):
            f_k, PSD_k[baseline] = signal.welch(opdc_kalman_pizeo_opd[:,baseline],
                    fs=(1./numpy.nanmean(numpy.diff(self.time))), detrend='linear',
                    nperseg=1024, scaling='spectrum')

            PSD_k[baseline] = numpy.sqrt(PSD_k[baseline])*1000.0
        self.PSD_k = PSD_k
        self.f_k = f_k

    def getPlateScales(self):
        retval = {}
        for i in self.AcqCamDat.SCALE.data.keys():
            retval[i] = numpy.mean(self.AcqCamDat.SCALE.data[i][numpy.isfinite(self.AcqCamDat.SCALE.data[i])])
        return retval

    def findVibrationPeaks(self, ax=None):
        self.flattenPSDs()
        self.peaks = {}
        for baseline in range(0,6):
            self.peaks[baseline] = {}
            peaks = scipy.signal.find_peaks_cwt(self.flattenedPSD_k[baseline],
                    numpy.arange(5,10))
            diffs = self.f_k[numpy.diff(peaks)]
            self.peaks[baseline]['freqs'] = []
            self.peaks[baseline]['power'] = []
            for i in range(len(diffs)):
                if diffs[i] > 3:
                    window = range(numpy.max([0, peaks[i]-3]), numpy.min([len(self.f_k), peaks[i]+3]))
                    self.peaks[baseline]['freqs'].append(self.f_k[peaks[i]])
                    self.peaks[baseline]['power'].append(scipy.integrate.trapz(self.PSD_k[baseline][window],
                                                 x=self.f_k[window]))

            self.peaks[baseline]['freqs'] = numpy.array(self.peaks[baseline]['freqs'])
            self.peaks[baseline]['power'] = numpy.array(self.peaks[baseline]['power'])
            
            if ax != None:
                #ax.clear()
                ax.set_xscale('linear')
                ax.set_yscale('linear')
                ax.plot(self.f_k, self.PSD_k[baseline])
                ax.figure.show()
                raw_input()

        return self.peaks

    def flattenPSDs(self):
        window = self.f_k > 15.0
        self.flattenedPSD_k = []
        for baseline in range(0,6):
            A = numpy.log10(self.PSD_k[baseline][window])
            model = numpy.array([numpy.ones(self.f_k[window].shape), numpy.log10(self.f_k[window])])
            im = numpy.linalg.pinv(model)
            SpectralSlope = im.T.dot(A.T)

            fullModel = numpy.array([numpy.ones(self.f_k.shape), numpy.log10(self.f_k)])
            self.flattenedPSD_k.append(self.PSD_k[baseline]/10.0**SpectralSlope.dot(fullModel))

class GRAVITY_Dual_Sci_P2VM( GRAVITY_Dual_P2VM ):
    def __init__(self, fileBase = '', startTime=0.0, CIAO_Data=None,
            processAcqCamData=False):
        self.filename = fileBase+'_dualscip2vmred.fits'
        self.startTime = startTime
        HDU = pyfits.open(self.filename)
        indices = []
        extnames = []
        for i, h in enumerate(HDU):
            extnames.append(h.name)
            indices.append(i)
        self.indices = numpy.array(indices)
        self.extnames = numpy.array(extnames)
        self.opdc = HDU['OPDC'].data
        self.opdc_kalman_piezo = self.opdc.field('KALMAN_PIEZO')  # Time, telescope
        self.opdc_kalman_opd = self.opdc.field('KALMAN_OPD')
        self.time = self.opdc.field('TIME')/1e6
        header = HDU["PRIMARY"].header
        self.polarization = header.get('ESO INS POLA MODE')
        SC = None
        FT = None
        if self.polarization == 'SPLIT':
            for oi in self.indices[self.extnames=='OI_FLUX']:
                if ('TOTALFLUX_FT' in HDU[oi].data.names):
                    BCI_Time = HDU[oi].data.field("TIME")
                    FDDL = HDU[oi].data.field("FDDL")
                    FT_POS = HDU[oi].data.field("FT_POS")
                    SC_POS = HDU[oi].data.field("SC_POS")
                    if SC == None:
                        SC = HDU[oi].data.field("TOTALFLUX_SC")
                        FT = HDU[oi].data.field("TOTALFLUX_FT")
                    else:
                        SC = SC + HDU[oi].data.field("TOTALFLUX_SC")
                        FT = FT + HDU[oi].data.field("TOTALFLUX_FT")
        elif self.polarization == 'COMBINED':
            for oi in self.indices[self.extnames=='OI_FLUX']:
                if ('TOTALFLUX_FT' in HDU[oi].data.names):
                    BCI_Time = HDU[oi].data.field("TIME")
                    FDDL = HDU[oi].data.field("FDDL")
                    FT_POS = HDU[oi].data.field("FT_POS")
                    SC_POS = HDU[oi].data.field("SC_POS")
                    SC = HDU[oi].data.field("TOTALFLUX_SC")
                    FT = HDU[oi].data.field("TOTALFLUX_FT")
        self.TOTALFLUX_SC = SC
        self.TOTALFLUX_FT = FT
        self.BCI_Time = BCI_Time
        self.FDDL = FDDL
        self.FT_POS = FT_POS
        self.SC_POS = SC_POS
        self.SObjMag = header.get('ESO INS SOBJ MAG')
        self.SObjName = header.get('ESO INS SOBJ NAME')
        self.AOInUse = header.get('ESO COU AO SYSTEM')
        self.AO_Mag = header.get('ESO COU GUID MAG')
        self.AO_Wave = header.get('ESO COU GUID WAVELEN')
        self.Kal_Gain = header.get('ESO FT KAL GAIN')
        self.Kal_Mode = header.get('ESO FT KAL MODE')
        self.Derot1 = header.get('ESO INS DROT1 ENC')
        self.Derot2 = header.get('ESO INS DROT2 ENC')
        self.Derot3 = header.get('ESO INS DROT3 ENC')
        self.Derot4 = header.get('ESO INS DROT4 ENC')
        self.Derot5 = header.get('ESO INS DROT5 ENC')
        self.Derot6 = header.get('ESO INS DROT6 ENC')
        self.Derot7 = header.get('ESO INS DROT7 ENC')
        self.Derot8 = header.get('ESO INS DROT8 ENC')
        self.SOBJ_OFFX = header.get('ESO INS SOBJ OFFX')
        self.SOBJ_OFFY = header.get('ESO INS SOBJ OFFY')
        self.MET_SOBJ_DRA = {}
        self.MET_SOBJ_DDEC = {}
        self.MET_SOBJ_DRA[0] = header.get('ESO QC MET SOBJ DRA1')
        self.MET_SOBJ_DDEC[0] = header.get('ESO QC MET SOBJ DDEC1')
        self.MET_SOBJ_DRA[1] = header.get('ESO QC MET SOBJ DRA2')
        self.MET_SOBJ_DDEC[1] = header.get('ESO QC MET SOBJ DDEC2')
        self.MET_SOBJ_DRA[2] = header.get('ESO QC MET SOBJ DRA3')
        self.MET_SOBJ_DDEC[2] = header.get('ESO QC MET SOBJ DDEC3')
        self.MET_SOBJ_DRA[3] = header.get('ESO QC MET SOBJ DRA4')
        self.MET_SOBJ_DDEC[3] = header.get('ESO QC MET SOBJ DDEC4')
        self.NORTH_ANGLE = {}
        self.NORTH_ANGLE[0] = header.get('ESO QC ACQ FIELD1 NORTH_ANGLE')
        self.NORTH_ANGLE[1] = header.get('ESO QC ACQ FIELD2 NORTH_ANGLE')
        self.NORTH_ANGLE[2] = header.get('ESO QC ACQ FIELD3 NORTH_ANGLE')
        self.NORTH_ANGLE[3] = header.get('ESO QC ACQ FIELD4 NORTH_ANGLE')
        self.Strehl = {}
        self.Strehl[0] = header.get('ESO QC ACQ FIELD1 STREHL')
        self.Strehl[1] = header.get('ESO QC ACQ FIELD2 STREHL')
        self.Strehl[2] = header.get('ESO QC ACQ FIELD3 STREHL')
        self.Strehl[3] = header.get('ESO QC ACQ FIELD4 STREHL')
        self.SObjSwap = 'YES' in header.get('ESO INS SOBJ SWAP')
        self.FTName = header.get('ESO FT ROBJ NAME')
        self.FTMag = header.get('ESO FT ROBJ MAG')
        self.AcqCamDIT = header.get('ESO DET1 SEQ1 DIT')
        self.ISS_Config = header.get('ESO ISS CONF T1NAME')[0]
        self.FAFT1_CURX = header.get('ESO INS FAFT1 CURX')
        self.FAFT1_CURY = header.get('ESO INS FAFT1 CURY')
        self.FAFT2_CURX = header.get('ESO INS FAFT2 CURX')
        self.FAFT2_CURY = header.get('ESO INS FAFT2 CURY')
        self.FAFT3_CURX = header.get('ESO INS FAFT3 CURX')
        self.FAFT3_CURY = header.get('ESO INS FAFT3 CURY')
        self.FAFT4_CURX = header.get('ESO INS FAFT4 CURX')
        self.FAFT4_CURY = header.get('ESO INS FAFT4 CURY')
        self.FAFT5_CURX = header.get('ESO INS FAFT5 CURX')
        self.FAFT5_CURY = header.get('ESO INS FAFT5 CURY')
        self.FAFT6_CURX = header.get('ESO INS FAFT6 CURX')
        self.FAFT6_CURY = header.get('ESO INS FAFT6 CURY')
        self.FAFT7_CURX = header.get('ESO INS FAFT7 CURX')
        self.FAFT7_CURY = header.get('ESO INS FAFT7 CURY')
        self.FAFT8_CURX = header.get('ESO INS FAFT8 CURX')
        self.FAFT8_CURY = header.get('ESO INS FAFT8 CURY')
        if processAcqCamData:
            for a in self.indices[self.extnames == 'IMAGING_DATA_ACQ']:
                if HDU[a].data.shape[0] > 1:
                    break
            imagingData = (header,HDU[a].data)
        else:
            imagingData = None
        #try:
        self.AcqCamDat = AcqCamData(fiberData=HDU[self.indices[self.extnames=="OI_VIS_ACQ"][0]].data,
                fluxData={"TOTALFLUX_FT":self.TOTALFLUX_FT,"TOTALFLUX_SC":self.TOTALFLUX_SC,"TIME":self.BCI_Time,"FDDL":self.FDDL,
                    "FT_POS": self.FT_POS, "SC_POS":self.SC_POS},
                FTMag=self.FTMag, 
                SCMag=self.SObjMag,
                AcqDit=self.AcqCamDIT,
                CIAO_Data=CIAO_Data,
                imagingData=imagingData)
        self.AcqCamDat.binData()
        self.AcqCamDat.calcMedian()
        #except:
        #    self.AcqCamDat = None

class GRAVITY_Dual_Cal_P2VM( GRAVITY_Dual_P2VM ):
    def __init__(self, fileBase = '', startTime=0.0, CIAO_Data=None,
            processAcqCamData=False):
        self.filename = fileBase+'_dualcalp2vmred.fits'
        self.startTime = startTime
        HDU = pyfits.open(self.filename)
        indices = []
        extnames = []
        for i, h in enumerate(HDU):
            extnames.append(h.name)
            indices.append(i)
        self.indices = numpy.array(indices)
        self.extnames = numpy.array(extnames)
        self.opdc = HDU['OPDC'].data
        self.opdc_kalman_piezo = self.opdc.field('KALMAN_PIEZO')  # Time, telescope
        self.opdc_kalman_opd = self.opdc.field('KALMAN_OPD')
        self.time = self.opdc.field('TIME')/1e6
        header = HDU["PRIMARY"].header
        self.polarization = header.get('ESO INS POLA MODE')
        SC = None
        FT = None
        if self.polarization == 'SPLIT':
            for oi in self.indices[self.extnames=='OI_FLUX']:
                if ('TOTALFLUX_FT' in HDU[oi].data.names):
                    BCI_Time = HDU[oi].data.field("TIME")
                    FDDL = HDU[oi].data.field("FDDL")
                    FT_POS = HDU[oi].data.field("FT_POS")
                    SC_POS = HDU[oi].data.field("SC_POS")
                    if SC == None:
                        SC = HDU[oi].data.field("TOTALFLUX_SC")
                        FT = HDU[oi].data.field("TOTALFLUX_FT")
                    else:
                        SC = SC + HDU[oi].data.field("TOTALFLUX_SC")
                        FT = FT + HDU[oi].data.field("TOTALFLUX_FT")
        elif self.polarization == 'COMBINED':
            for oi in self.indices[self.extnames=='OI_FLUX']:
                if ('TOTALFLUX_FT' in HDU[oi].data.names):
                    BCI_Time = HDU[oi].data.field("TIME")
                    FDDL = HDU[oi].data.field("FDDL")
                    FT_POS = HDU[oi].data.field("FT_POS")
                    SC_POS = HDU[oi].data.field("SC_POS")
                    SC = HDU[oi].data.field("TOTALFLUX_SC")
                    FT = HDU[oi].data.field("TOTALFLUX_FT")
        self.TOTALFLUX_SC = SC
        self.TOTALFLUX_FT = FT
        self.BCI_Time = BCI_Time
        self.FDDL = FDDL
        self.FT_POS = FT_POS
        self.SC_POS = SC_POS
        self.SObjMag = header.get('ESO INS SOBJ MAG')
        self.SObjName = header.get('ESO INS SOBJ NAME')
        self.AOInUse = header.get('ESO COU AO SYSTEM')
        self.AO_Mag = header.get('ESO COU GUID MAG')
        self.AO_Wave = header.get('ESO COU GUID WAVELEN')
        self.Kal_Gain = header.get('ESO FT KAL GAIN')
        self.Kal_Mode = header.get('ESO FT KAL MODE')
        self.Derot1 = header.get('ESO INS DROT1 ENC')
        self.Derot2 = header.get('ESO INS DROT2 ENC')
        self.Derot3 = header.get('ESO INS DROT3 ENC')
        self.Derot4 = header.get('ESO INS DROT4 ENC')
        self.Derot5 = header.get('ESO INS DROT5 ENC')
        self.Derot6 = header.get('ESO INS DROT6 ENC')
        self.Derot7 = header.get('ESO INS DROT7 ENC')
        self.Derot8 = header.get('ESO INS DROT8 ENC')
        self.SOBJ_OFFX = header.get('ESO INS SOBJ OFFX')
        self.SOBJ_OFFY = header.get('ESO INS SOBJ OFFY')
        self.SObjSwap = 'YES' in header.get('ESO INS SOBJ SWAP')
        self.FTName = header.get('ESO FT ROBJ NAME')
        self.FTMag = header.get('ESO FT ROBJ MAG')
        self.AcqCamDIT = header.get('ESO DET1 SEQ1 DIT')
        self.ISS_Config = 'CAL'
        self.FAFT1_CURX = header.get('ESO INS FAFT1 CURX')
        self.FAFT1_CURY = header.get('ESO INS FAFT1 CURY')
        self.FAFT2_CURX = header.get('ESO INS FAFT2 CURX')
        self.FAFT2_CURY = header.get('ESO INS FAFT2 CURY')
        self.FAFT3_CURX = header.get('ESO INS FAFT3 CURX')
        self.FAFT3_CURY = header.get('ESO INS FAFT3 CURY')
        self.FAFT4_CURX = header.get('ESO INS FAFT4 CURX')
        self.FAFT4_CURY = header.get('ESO INS FAFT4 CURY')
        self.FAFT5_CURX = header.get('ESO INS FAFT5 CURX')
        self.FAFT5_CURY = header.get('ESO INS FAFT5 CURY')
        self.FAFT6_CURX = header.get('ESO INS FAFT6 CURX')
        self.FAFT6_CURY = header.get('ESO INS FAFT6 CURY')
        self.FAFT7_CURX = header.get('ESO INS FAFT7 CURX')
        self.FAFT7_CURY = header.get('ESO INS FAFT7 CURY')
        self.FAFT8_CURX = header.get('ESO INS FAFT8 CURX')
        self.FAFT8_CURY = header.get('ESO INS FAFT8 CURY')
        if processAcqCamData:
            for a in self.indices[self.extnames == 'IMAGING_DATA_ACQ']:
                if HDU[a].data.shape[0] > 1:
                    break
            imagingData = (header,HDU[a].data)
        else:
            imagingData = None
        try:
            self.AcqCamDat = AcqCamData(fiberData=HDU["OI_VIS_ACQ"].data,
                fluxData={"TOTALFLUX_FT":self.TOTALFLUX_FT,"TOTALFLUX_SC":self.TOTALFLUX_SC,"TIME":self.BCI_Time,"FDDL":self.FDDL,
                    "FT_POS": self.FT_POS, "SC_POS":self.SC_POS},
                FTMag=self.FTMag, 
                SCMag=self.SObjMag,
                AcqDit=self.AcqCamDIT,
                CIAO_Data=CIAO_Data,
                imagingData=imagingData)
            self.AcqCamDat.binData()
            self.AcqCamDat.calcMedian()
        except:
            self.AcqCamDat = None

class GRAVITY_Single_P2VM ( object ):
    def __init__(self, fileBase = '', startTime=0.0, CIAO_Data=None,
            processAcqCamData=False):
        return self

    def getPupilMotion(self):
        self.TIME = (self.data.field('TIME')[::4]*1e-6)/86400.0+self.startTime.mjd
        self.PUPIL_NSPOT = self.data.field('PUPIL_NSPOT')[::4]
        self.PUPIL_X = {}
        self.PUPIL_Y = {}
        self.PUPIL_Z = {}
        self.PUPIL_R = {}
        self.PUPIL_U = {}
        self.PUPIL_V = {}
        self.PUPIL_W = {}
        
        for UT in [1, 2, 3, 4]:
            self.PUPIL_X[UT] = self.data.field('PUPIL_X')[UT-1::4]
            self.PUPIL_Y[UT] = self.data.field('PUPIL_Y')[UT-1::4]
            self.PUPIL_Z[UT] = self.data.field('PUPIL_Z')[UT-1::4]
            self.PUPIL_R[UT] = self.data.field('PUPIL_R')[UT-1::4]
            self.PUPIL_U[UT] = self.data.field('PUPIL_U')[UT-1::4]
            self.PUPIL_V[UT] = self.data.field('PUPIL_V')[UT-1::4]
            self.PUPIL_W[UT] = self.data.field('PUPIL_W')[UT-1::4]

    def computeOPDPeriodograms(self):
        M_matrix = numpy.array([-1.,1.,0.0,0.0,-1.,0.0,1.,0.0,-1.,0.0,0.0,1.,0.0,-1.,1.,0.0,0.0,-1.,0.0,1.,0.0,0.0,-1.,1.]);
        M_matrix = M_matrix.reshape((6, 4))

        opdc_kalman_pizeo_opd = numpy.dot(M_matrix,
                self.opdc_kalman_piezo.T).T
        PSD_k = [1, 2, 3, 4, 5, 6]
        for baseline in range(0, 6):
            f_k, PSD_k[baseline] = signal.welch(opdc_kalman_pizeo_opd[:,baseline],
                    fs=(1./numpy.nanmean(numpy.diff(self.time))), detrend='linear',
                    nperseg=1024, scaling='spectrum')

            PSD_k[baseline] = numpy.sqrt(PSD_k[baseline])*1000.0
        self.PSD_k = PSD_k
        self.f_k = f_k

    def findVibrationPeaks(self, ax=None):
        self.flattenPSDs()
        self.peaks = {}
        for baseline in range(0,6):
            self.peaks[baseline] = {}
            peaks = scipy.signal.find_peaks_cwt(self.flattenedPSD_k[baseline],
                    numpy.arange(5,10))
            diffs = self.f_k[numpy.diff(peaks)]
            self.peaks[baseline]['freqs'] = []
            self.peaks[baseline]['power'] = []
            for i in range(len(diffs)):
                if diffs[i] > 3:
                    window = range(numpy.max([0, peaks[i]-3]), numpy.min([len(self.f_k), peaks[i]+3]))
                    self.peaks[baseline]['freqs'].append(self.f_k[peaks[i]])
                    self.peaks[baseline]['power'].append(scipy.integrate.trapz(self.PSD_k[baseline][window],
                                                 x=self.f_k[window]))

            self.peaks[baseline]['freqs'] = numpy.array(self.peaks[baseline]['freqs'])
            self.peaks[baseline]['power'] = numpy.array(self.peaks[baseline]['power'])
            
            if ax != None:
                #ax.clear()
                ax.set_xscale('linear')
                ax.set_yscale('linear')
                ax.plot(self.f_k, self.PSD_k[baseline])
                ax.figure.show()
                raw_input()

        return self.peaks

    def getPlateScales(self):
        retval = {}
        for i in self.AcqCamDat.SCALE.data.keys():
            retval[i] = numpy.mean(self.AcqCamDat.SCALE.data[i][numpy.isfinite(self.AcqCamDat.SCALE.data[i])])
        return retval



    def flattenPSDs(self):
        window = self.f_k > 15.0
        self.flattenedPSD_k = []
        for baseline in range(0,6):
            A = numpy.log10(self.PSD_k[baseline][window])
            model = numpy.array([numpy.ones(self.f_k[window].shape), numpy.log10(self.f_k[window])])
            im = numpy.linalg.pinv(model)
            SpectralSlope = im.T.dot(A.T)

            fullModel = numpy.array([numpy.ones(self.f_k.shape), numpy.log10(self.f_k)])
            self.flattenedPSD_k.append(self.PSD_k[baseline]/10.0**SpectralSlope.dot(fullModel))


class GRAVITY_Single_Sci_P2VM( GRAVITY_Single_P2VM ):
    def __init__(self, fileBase = '', startTime=0.0, CIAO_Data=None,
            processAcqCamData=False):
        self.filename = fileBase+'_singlescip2vmred.fits'
        self.startTime = startTime
        HDU = pyfits.open(self.filename)
        indices = []
        extnames = []
        for i, h in enumerate(HDU):
            extnames.append(h.name)
            indices.append(i)
        self.indices = numpy.array(indices)
        self.extnames = numpy.array(extnames)
        self.opdc = HDU['OPDC'].data
        self.opdc_kalman_piezo = self.opdc.field('KALMAN_PIEZO')  # Time, telescope
        self.opdc_kalman_opd = self.opdc.field('KALMAN_OPD')
        self.time = self.opdc.field('TIME')/1e6
        header = HDU["PRIMARY"].header
        self.polarization = header.get('ESO INS POLA MODE')
        SC = None
        FT = None
        if self.polarization == 'SPLIT':
            for oi in self.indices[self.extnames=='OI_FLUX']:
                if ('TOTALFLUX_FT' in HDU[oi].data.names):
                    BCI_Time = HDU[oi].data.field("TIME")
                    FDDL = HDU[oi].data.field("FDDL")
                    FT_POS = HDU[oi].data.field("FT_POS")
                    SC_POS = HDU[oi].data.field("SC_POS")
                    if SC == None:
                        SC = HDU[oi].data.field("TOTALFLUX_SC")
                        FT = HDU[oi].data.field("TOTALFLUX_FT")
                    else:
                        SC = SC + HDU[oi].data.field("TOTALFLUX_SC")
                        FT = FT + HDU[oi].data.field("TOTALFLUX_FT")
        elif self.polarization == 'COMBINED':
            for oi in self.indices[self.extnames=='OI_FLUX']:
                if ('TOTALFLUX_FT' in HDU[oi].data.names):
                    BCI_Time = HDU[oi].data.field("TIME")
                    FDDL = HDU[oi].data.field("FDDL")
                    FT_POS = HDU[oi].data.field("FT_POS")
                    SC_POS = HDU[oi].data.field("SC_POS")
                    SC = HDU[oi].data.field("TOTALFLUX_SC")
                    FT = HDU[oi].data.field("TOTALFLUX_FT")
        self.TOTALFLUX_SC = SC
        self.TOTALFLUX_FT = FT
        self.BCI_Time = BCI_Time
        self.FDDL = FDDL
        self.FT_POS = FT_POS
        self.SC_POS = SC_POS
        self.SObjMag = header.get('ESO INS SOBJ MAG')
        self.SObjName = header.get('ESO INS SOBJ NAME')
        self.AOInUse = header.get('ESO COU AO SYSTEM')
        self.AO_Mag = header.get('ESO COU GUID MAG')
        self.AO_Wave = header.get('ESO COU GUID WAVELEN')
        self.Kal_Gain = header.get('ESO FT KAL GAIN')
        self.Kal_Mode = header.get('ESO FT KAL MODE')
        self.Derot1 = header.get('ESO INS DROT1 ENC')
        self.Derot2 = header.get('ESO INS DROT2 ENC')
        self.Derot3 = header.get('ESO INS DROT3 ENC')
        self.Derot4 = header.get('ESO INS DROT4 ENC')
        self.Derot5 = header.get('ESO INS DROT5 ENC')
        self.Derot6 = header.get('ESO INS DROT6 ENC')
        self.Derot7 = header.get('ESO INS DROT7 ENC')
        self.Derot8 = header.get('ESO INS DROT8 ENC')
        self.SOBJ_OFFX = header.get('ESO INS SOBJ OFFX')
        self.SOBJ_OFFY = header.get('ESO INS SOBJ OFFY')
        swap = header.get('ESO INS SOBJ SWAP')
        self.ISS_Config = header.get('ESO ISS CONF T1NAME')[0]
        self.FAFT1_CURX = header.get('ESO INS FAFT1 CURX')
        self.FAFT1_CURY = header.get('ESO INS FAFT1 CURY')
        self.FAFT2_CURX = header.get('ESO INS FAFT2 CURX')
        self.FAFT2_CURY = header.get('ESO INS FAFT2 CURY')
        self.FAFT3_CURX = header.get('ESO INS FAFT3 CURX')
        self.FAFT3_CURY = header.get('ESO INS FAFT3 CURY')
        self.FAFT4_CURX = header.get('ESO INS FAFT4 CURX')
        self.FAFT4_CURY = header.get('ESO INS FAFT4 CURY')
        self.FAFT5_CURX = header.get('ESO INS FAFT5 CURX')
        self.FAFT5_CURY = header.get('ESO INS FAFT5 CURY')
        self.FAFT6_CURX = header.get('ESO INS FAFT6 CURX')
        self.FAFT6_CURY = header.get('ESO INS FAFT6 CURY')
        self.FAFT7_CURX = header.get('ESO INS FAFT7 CURX')
        self.FAFT7_CURY = header.get('ESO INS FAFT7 CURY')
        self.FAFT8_CURX = header.get('ESO INS FAFT8 CURX')
        self.FAFT8_CURY = header.get('ESO INS FAFT8 CURY')
        if swap != None:
            self.SObjSwap = 'YES' in header.get('ESO INS SOBJ SWAP')
        else:
            self.SObjSwap = False
        self.FTName = header.get('ESO FT ROBJ NAME')
        self.FTMag = header.get('ESO FT ROBJ MAG')
        self.AcqCamDIT = header.get('ESO DET1 SEQ1 DIT')
        if processAcqCamData:
            for a in self.indices[self.extnames == 'IMAGING_DATA_ACQ']:
                if HDU[a].data.shape[0] > 1:
                    break
            imagingData = (header,HDU[a].data)
        else:
            imagingData = None
        try:
            self.AcqCamDat = AcqCamData(fiberData=HDU["OI_VIS_ACQ"].data,
                fluxData={"TOTALFLUX_FT":self.TOTALFLUX_FT,"TOTALFLUX_SC":self.TOTALFLUX_SC,"TIME":self.BCI_Time,"FDDL":self.FDDL,
                    "FT_POS": self.FT_POS, "SC_POS":self.SC_POS},
                FTMag=self.FTMag, 
                SCMag=self.SObjMag,
                AcqDit=self.AcqCamDIT,
                CIAO_Data=CIAO_Data,
                imagingData=imagingData)
            self.AcqCamDat.binData()
            self.AcqCamDat.calcMedian()
            #self.getPupilMotion()
        except:
            self.AcqCamDat = None

class GRAVITY_Single_Cal_P2VM( GRAVITY_Single_P2VM ):
    def __init__(self, fileBase = '', startTime=0.0, CIAO_Data=None,
            processAcqCamData=False):
        self.filename = fileBase+'_singlecalp2vmred.fits'
        self.startTime = startTime
        HDU = pyfits.open(self.filename)
        indices = []
        extnames = []
        for i, h in enumerate(HDU):
            extnames.append(h.name)
            indices.append(i)
        self.indices = numpy.array(indices)
        self.extnames = numpy.array(extnames)
        self.opdc = HDU['OPDC'].data
        self.opdc_kalman_piezo = self.opdc.field('KALMAN_PIEZO')  # Time, telescope
        self.opdc_kalman_opd = self.opdc.field('KALMAN_OPD')
        self.time = self.opdc.field('TIME')/1e6
        header = HDU["PRIMARY"].header
        self.polarization = header.get('ESO INS POLA MODE')
        SC = None
        FT = None
        if self.polarization == 'SPLIT':
            for oi in self.indices[self.extnames=='OI_FLUX']:
                if ('TOTALFLUX_FT' in HDU[oi].data.names):
                    BCI_Time = HDU[oi].data.field("TIME")
                    FDDL = HDU[oi].data.field("FDDL")
                    FT_POS = HDU[oi].data.field("FT_POS")
                    SC_POS = HDU[oi].data.field("SC_POS")
                    if SC == None:
                        SC = HDU[oi].data.field("TOTALFLUX_SC")
                        FT = HDU[oi].data.field("TOTALFLUX_FT")
                    else:
                        SC = SC + HDU[oi].data.field("TOTALFLUX_SC")
                        FT = FT + HDU[oi].data.field("TOTALFLUX_FT")
        elif self.polarization == 'COMBINED':
            for oi in self.indices[self.extnames=='OI_FLUX']:
                if ('TOTALFLUX_FT' in HDU[oi].data.names):
                    BCI_Time = HDU[oi].data.field("TIME")
                    FDDL = HDU[oi].data.field("FDDL")
                    FT_POS = HDU[oi].data.field("FT_POS")
                    SC_POS = HDU[oi].data.field("SC_POS")
                    SC = HDU[oi].data.field("TOTALFLUX_SC")
                    FT = HDU[oi].data.field("TOTALFLUX_FT")
        self.TOTALFLUX_SC = SC
        self.TOTALFLUX_FT = FT
        self.BCI_Time = BCI_Time
        self.FDDL = FDDL
        self.FT_POS = FT_POS
        self.SC_POS = SC_POS
        self.SObjMag = header.get('ESO INS SOBJ MAG')
        self.SObjName = header.get('ESO INS SOBJ NAME')
        self.AOInUse = header.get('ESO COU AO SYSTEM')
        self.AO_Mag = header.get('ESO COU GUID MAG')
        self.AO_Wave = header.get('ESO COU GUID WAVELEN')
        self.Kal_Gain = header.get('ESO FT KAL GAIN')
        self.Kal_Mode = header.get('ESO FT KAL MODE')
        self.Derot1 = header.get('ESO INS DROT1 ENC')
        self.Derot2 = header.get('ESO INS DROT2 ENC')
        self.Derot3 = header.get('ESO INS DROT3 ENC')
        self.Derot4 = header.get('ESO INS DROT4 ENC')
        self.Derot5 = header.get('ESO INS DROT5 ENC')
        self.Derot6 = header.get('ESO INS DROT6 ENC')
        self.Derot7 = header.get('ESO INS DROT7 ENC')
        self.Derot8 = header.get('ESO INS DROT8 ENC')
        self.SOBJ_OFFX = header.get('ESO INS SOBJ OFFX')
        self.SOBJ_OFFY = header.get('ESO INS SOBJ OFFY')
        swap = header.get('ESO INS SOBJ SWAP')
        self.ISS_Config = 'CAL'
        self.FAFT1_CURX = header.get('ESO INS FAFT1 CURX')
        self.FAFT1_CURY = header.get('ESO INS FAFT1 CURY')
        self.FAFT2_CURX = header.get('ESO INS FAFT2 CURX')
        self.FAFT2_CURY = header.get('ESO INS FAFT2 CURY')
        self.FAFT3_CURX = header.get('ESO INS FAFT3 CURX')
        self.FAFT3_CURY = header.get('ESO INS FAFT3 CURY')
        self.FAFT4_CURX = header.get('ESO INS FAFT4 CURX')
        self.FAFT4_CURY = header.get('ESO INS FAFT4 CURY')
        self.FAFT5_CURX = header.get('ESO INS FAFT5 CURX')
        self.FAFT5_CURY = header.get('ESO INS FAFT5 CURY')
        self.FAFT6_CURX = header.get('ESO INS FAFT6 CURX')
        self.FAFT6_CURY = header.get('ESO INS FAFT6 CURY')
        self.FAFT7_CURX = header.get('ESO INS FAFT7 CURX')
        self.FAFT7_CURY = header.get('ESO INS FAFT7 CURY')
        self.FAFT8_CURX = header.get('ESO INS FAFT8 CURX')
        self.FAFT8_CURY = header.get('ESO INS FAFT8 CURY')
        if swap != None:
            self.SObjSwap = 'YES' in header.get('ESO INS SOBJ SWAP')
        else:
            self.SObjSwap = False
        self.FTName = header.get('ESO FT ROBJ NAME')
        self.FTMag = header.get('ESO FT ROBJ MAG')
        self.AcqCamDIT = header.get('ESO DET1 SEQ1 DIT')
        if processAcqCamData:
            for a in self.indices[self.extnames == 'IMAGING_DATA_ACQ']:
                if HDU[a].data.shape[0] > 1:
                    break
            imagingData = (header,HDU[a].data)
        else:
            imagingData = None
        try:
            self.AcqCamDat = AcqCamData(fiberData=HDU["OI_VIS_ACQ"].data,
                fluxData={"TOTALFLUX_FT":self.TOTALFLUX_FT,"TOTALFLUX_SC":self.TOTALFLUX_SC,"TIME":self.BCI_Time,"FDDL":self.FDDL,
                    "FT_POS": self.FT_POS, "SC_POS":self.SC_POS},
                FTMag=self.FTMag, 
                SCMag=self.SObjMag,
                AcqDit=self.AcqCamDIT,
                CIAO_Data=CIAO_Data,
                imagingData=imagingData)
            self.AcqCamDat.binData()
            self.AcqCamDat.calcMedian()
            #self.getPupilMotion()
        except:
            self.AcqCamDat = None


class GRAVITY_Dual_SciVis( object ):
    def __init__(self, fileBase = ''):
        self.filename = fileBase+'_dualscivis.fits'
        try:
            self.header = pyfits.getheader(self.filename)
            self.data = pyfits.getdata(self.filename)
        except:
            pass

class GRAVITY_AstroReduced( object ):
    def __init__(self, fileBase = ''):
        self.filename = fileBase+'_astroreduced.fits'
        self.masterheader = pyfits.getheader(self.filename, ext=0)
        self.RA = self.masterheader.get('RA')
        self.DEC = self.masterheader.get('DEC')
        HDUList = pyfits.open(self.filename)
        self.data = []
        self.headers = []
        for HDU in HDUList:
            if HDU.name == 'OI_FLUX':
                self.data.append(HDU.data)
                self.headers.append(HDU.header)
        HDUList.close()

class GRAVITY_Data( object ):
    def __init__(self, fileBase='', UTS=[1,2,3,4], sqlCursor=None, CIAO_Data =
            None, processAcqCamData=False, datadir=None):
        if datadir != None:
            self.datadir = datadir
        else:
            self.datadir = os.environ.get('GRAVITY_DATA')
        if self.datadir in fileBase:
            self.fileBase = str(fileBase[len(self.datadir):])
        else:
            self.fileBase = fileBase
        self.sqlCursor = sqlCursor
        self.UTS = UTS
        #self.Dual_SciVis = GRAVITY_Dual_SciVis(fileBase=self.datadir+self.fileBase)
        self.AstroReduced = GRAVITY_AstroReduced(fileBase=self.datadir+self.fileBase)
        #self.startTime = time.mktime(datetime.strptime(self.AstroReduced.masterheader.get('ESO PCR ACQ START'), '%Y-%m-%dT%H:%M:%S.%f').timetuple())
        self.startTime = aptime.Time(self.AstroReduced.masterheader.get('ESO PCR ACQ START'))
        self.AcqCamData = None
        self.SingleSciP2VM = None
        self.DualSciP2VM = None
        self.SingleCalP2VM = None
        self.DualCalP2VM = None
        if os.path.exists(self.datadir+self.fileBase+'_singlescip2vmred.fits'):
            self.SingleSciP2VM = GRAVITY_Single_Sci_P2VM(fileBase=self.datadir+self.fileBase,
                startTime=self.startTime, CIAO_Data=CIAO_Data,
                processAcqCamData=processAcqCamData)
            self.AcqCamData = self.SingleSciP2VM.AcqCamDat
        elif os.path.exists(self.datadir+self.fileBase+'_dualscip2vmred.fits'):
            self.DualSciP2VM = GRAVITY_Dual_Sci_P2VM(fileBase=self.datadir+self.fileBase,
                startTime=self.startTime, CIAO_Data=CIAO_Data,
                processAcqCamData=processAcqCamData)
            self.AcqCamData = self.DualSciP2VM.AcqCamDat
        elif os.path.exists(self.datadir+self.fileBase+'_singlecalp2vmred.fits'):
            self.SingleCalP2VM = GRAVITY_Single_Sci_P2VM(fileBase=self.datadir+self.fileBase,
                startTime=self.startTime, CIAO_Data=CIAO_Data,
                processAcqCamData=processAcqCamData)
            self.AcqCamData = self.SingleCalP2VM.AcqCamDat
        elif os.path.exists(self.datadir+self.fileBase+'_dualcalp2vmred.fits'):
            self.DualCalP2VM = GRAVITY_Dual_Cal_P2VM(fileBase=self.datadir+self.fileBase,
                startTime=self.startTime, CIAO_Data=CIAO_Data,
                processAcqCamData=processAcqCamData)
            self.AcqCamData = self.DualCalP2VM.AcqCamDat
        self.DITTimes = {}
        self.SC_Fluxes = {}
        self.FT_Fluxes = {}
        self.Strehl = {}
        self.Seeing = {}
        self.ASM_Seeing = {}
        for UT in self.UTS:
            self.DITTimes[UT] = []
            self.SC_Fluxes[UT] = []
            self.FT_Fluxes[UT] = []
            for d in self.AstroReduced.data:
                self.DITTimes[UT].append((d.field('TIME')[UT-1::4]*1e-6)/(24*3600.0)+self.startTime.mjd)
                self.SC_Fluxes[UT].append(d.field('TOTALFLUX_SC')[UT-1::4])
                self.FT_Fluxes[UT].append(d.field('TOTALFLUX_FT')[UT-1::4])

    def computeOPDPeriodograms(self):
        if self.SingleSciP2VM != None:
            self.SingleSciP2VM.computeOPDPeriodograms()
        if self.DualSciP2VM != None:
            self.DualSciP2VM.computeOPDPeriodograms()

    def findVibrationPeaks(self):
        if self.SingleSciP2VM != None:
            return self.SingleSciP2VM.findVibrationPeaks()
        if self.DualSciP2VM != None:
            return self.DualSciP2VM.findVibrationPeaks()

    def getDerotatorPositions(self):
        derot = []
        if self.SingleSciP2VM != None:
            derot.append(self.SingleSciP2VM.Derot1)
            derot.append(self.SingleSciP2VM.Derot2)
            derot.append(self.SingleSciP2VM.Derot3)
            derot.append(self.SingleSciP2VM.Derot4)
            derot.append(self.SingleSciP2VM.Derot5)
            derot.append(self.SingleSciP2VM.Derot6)
            derot.append(self.SingleSciP2VM.Derot7)
            derot.append(self.SingleSciP2VM.Derot8)
        if self.DualSciP2VM != None:
            derot.append(self.DualSciP2VM.Derot1)
            derot.append(self.DualSciP2VM.Derot2)
            derot.append(self.DualSciP2VM.Derot3)
            derot.append(self.DualSciP2VM.Derot4)
            derot.append(self.DualSciP2VM.Derot5)
            derot.append(self.DualSciP2VM.Derot6)
            derot.append(self.DualSciP2VM.Derot7)
            derot.append(self.DualSciP2VM.Derot8)

        return numpy.array(derot)


    def getSobj_Offsets(self):
        if self.SingleSciP2VM != None:
            return (self.SingleSciP2VM.SOBJ_OFFX, self.SingleSciP2VM.SOBJ_OFFY)
        if self.DualSciP2VM != None:
            return (self.DualSciP2VM.SOBJ_OFFX, self.DualSciP2VM.SOBJ_OFFY)

    def getArrayConfig(self):
        if self.SingleSciP2VM != None:
            return self.SingleSciP2VM.ISS_Config
        if self.DualSciP2VM != None:
            return self.DualSciP2VM.ISS_Config

    def getPlateScales(self):
        if self.SingleSciP2VM != None:
            return self.SingleSciP2VM.getPlateScales()
        if self.DualSciP2VM != None:
            return self.DualSciP2VM.getPlateScales()
        

    def computeACQCAMStrehl(self):
        print asdf
        #AcqCamImage = AcqCamImage(
        
    def addToDatabase(self):
        if self.sqlCursor == None:
            return
        if self.AcqCamData == False:
            return
        sqlCommand = "SELECT * FROM GRAVITY_OBS WHERE TIMESTAMP = %.7f" % self.startTime.mjd
        self.sqlCursor.execute(sqlCommand)
        result = self.sqlCursor.fetchall()
        if len(result) == 0:
            values = {}
            values["TIMESTAMP"] = self.startTime.mjd
            values["DIRECTORY"] = self.datadir+self.fileBase
            
            if self.DualSciP2VM != None:
                values["FTOBJ_NAME"] = self.DualSciP2VM.FTName
                values["FTMAG"] = self.DualSciP2VM.FTMag
                values["SOBJ_NAME"] = self.DualSciP2VM.SObjName
                values["SOBJMAG"] = self.DualSciP2VM.SObjMag
                values["AO_SYS"] = self.DualSciP2VM.AOInUse
                values["AO_MAG"] = self.DualSciP2VM.AO_Mag
                values["AO_WAVE"] = self.DualSciP2VM.AO_Wave
                values["KAL_GAIN"] = self.DualSciP2VM.Kal_Gain
                values["KAL_MODE"] = self.DualSciP2VM.Kal_Mode
                values["DEROT1"] = self.DualSciP2VM.Derot1
                values["DEROT2"] = self.DualSciP2VM.Derot2
                values["DEROT3"] = self.DualSciP2VM.Derot3
                values["DEROT4"] = self.DualSciP2VM.Derot4
                values["DEROT5"] = self.DualSciP2VM.Derot5
                values["DEROT6"] = self.DualSciP2VM.Derot6
                values["DEROT7"] = self.DualSciP2VM.Derot7
                values["DEROT8"] = self.DualSciP2VM.Derot8
            elif self.SingleSciP2VM != None:
                values["FTOBJ_NAME"] = self.SingleSciP2VM.FTName
                values["FTMAG"] = self.SingleSciP2VM.FTMag
                values["SOBJ_NAME"] = self.SingleSciP2VM.SObjName
                values["SOBJMAG"] = self.SingleSciP2VM.SObjMag
                values["AO_SYS"] = self.SingleSciP2VM.AOInUse
                values["AO_MAG"] = self.SingleSciP2VM.AO_Mag
                values["AO_WAVE"] = self.SingleSciP2VM.AO_Wave
                values["KAL_GAIN"] = self.SingleSciP2VM.Kal_Gain
                values["KAL_MODE"] = self.SingleSciP2VM.Kal_Mode
                values["DEROT1"] = self.SingleSciP2VM.Derot1
                values["DEROT2"] = self.SingleSciP2VM.Derot2
                values["DEROT3"] = self.SingleSciP2VM.Derot3
                values["DEROT4"] = self.SingleSciP2VM.Derot4
                values["DEROT5"] = self.SingleSciP2VM.Derot5
                values["DEROT6"] = self.SingleSciP2VM.Derot6
                values["DEROT7"] = self.SingleSciP2VM.Derot7
                values["DEROT8"] = self.SingleSciP2VM.Derot8
            else:
                print "Oops! Didn't work!"
                return


            format_str= """INSERT INTO GRAVITY_OBS (TIMESTAMP, DIRECTORY,
            AO_SYSTEM, AO_GUIDE_MAG, AO_GUIDE_WAVELENGTH, FTOBJ_NAME, FTMAG,
            SOBJ_NAME, SOBJMAG, KALMAN_GAIN, KALMAN_MODE, DEROT1, DEROT2,
            DEROT3, DEROT4, DEROT5, DEROT6, DEROT7, DEROT8)  
                 VALUES ('{timestamp}', '{directory}', '{ao_sys}', '{ao_mag}',
                 '{ao_wave}', '{ftname}', '{ftmag}', '{sobjname}', '{sobjmag}',
                 '{kal_gain}', '{kal_mode}', '{derot1}', '{derot2}', '{derot3}',
                 '{derot4}', '{derot5}', '{derot6}', '{derot7}', '{derot8}');"""
                
            sqlCommand=format_str.format(timestamp=values['TIMESTAMP'],directory=values['DIRECTORY'],
                    ao_sys=values["AO_SYS"], ao_mag=values['AO_MAG'],
                    ao_wave=values['AO_WAVE'], ftname=values['FTOBJ_NAME'], ftmag=values['FTMAG'],
                    sobjname=values['SOBJ_NAME'], sobjmag=values['SOBJMAG'],
                    kal_gain=values['KAL_GAIN'], kal_mode=values['KAL_MODE'],
                    derot1=values["DEROT1"], derot2=values["DEROT2"],
                    derot3=values["DEROT3"], derot4=values["DEROT4"],
                    derot5=values["DEROT5"], derot6=values["DEROT6"],
                    derot7=values["DEROT7"], derot8=values["DEROT8"])
            
            self.sqlCursor.execute(sqlCommand)


    def plotReducedSCFluxes(self, ax=None):
        ax.clear()
        for i in self.UTS:
            ax.plot(self.DITTimes[i], self.SC_Fluxes[i])
        ax.figure.show()

    def plotReducedFTFluxes(self, ax=None):
        ax.clear()
        for i in range(4):
            ax.plot(self.DITTimes[i], self.FT_Fluxes[i])
        ax.figure.show()

    def addAOData(self, UT=0, Strehl = None, Seeing = None, ASM_Seeing = None):
        if Strehl != None:
            self.Strehl[UT] = Strehl
        if Seeing != None:
            self.Seeing[UT] = Seeing
        if ASM_Seeing != None:
            self.ASM_Seeing[UT] = ASM_Seeing

    def __lt__(self, other):
        if hasattr(other, startTime):
            return self.startTime < other.startTime
        else:
            return self.startTime < other

    def __le__(self, other):
        if hasattr(other, startTime):
            return self.startTime <= other.startTime
        else:
            return self.startTime <= other

    def __gt__(self, other):
        if hasattr(other, startTime):
            return self.startTime > other.startTime
        else:
            return self.startTime > other

    def __ge__(self, other):
        if hasattr(other, startTime):
            return self.startTime >= other.startTime
        else:
            return self.startTime >= other

    def __eq__(self, other):
        if hasattr(other, startTime):
            return self.startTime == other.startTime
        else:
            return self.startTime == other

    def __ne__(self, other):
        if hasattr(other, startTime):
            return self.startTime != other.startTime
        else:
            return self.startTime != other


class Anamorph( object ):
    def __init__(self, directory='', index=0):
        self.datadir = os.environ.get('CIAO_DATA')+directory
        self.CDMS_BaseDir = os.environ.get('CDMS_BASEDIR')
        self.CDMS_ConfigDir = os.environ.get('CDMS_CONFIGDIR')
        self.index = index

    def loadData(self):
        self.header = pyfits.getheader(self.datadir+'/CIAO_LOOP_%04d.fits'%self.index)
        self.CIAO_ID = self.header.get("ESO OCS SYS ID")
        loopData = pyfits.getdata(self.datadir+'/CIAO_LOOP_%04d.fits'%self.index)
        self.Intensities = loopData.field('Intensities')
        self.Gradients = loopData.field('Gradients')
        self.HODM = loopData.field('HODM_Positions')
        self.TTM = loopData.field('ITTM_Positions')
        self.FrameCounter = loopData.field('FrameCounter')
        self.startTime = time.mktime(time.strptime(self.header.get('ESO TPL START'),
                                     '%Y-%m-%dT%H:%M:%S'))
        self.loopTime = loopData.field('Seconds')+loopData.field('USeconds')/100000.0
        self.loopTime -= self.loopTime[0]
        self.loopRate = self.header.get('ESO AOS LOOP RATE')
        self.controlGain = self.header.get('ESO AOS GLOBAL GAIN')
        pixelData = pyfits.getdata(self.datadir+'/CIAO_PIXELS_%04d.fits'%self.index)
        self.Pixels = pixelData.field('Pixels')
        self.S2M = pyfits.getdata(self.datadir+'/RecnOptimiser.S2M_%04d.fits'%self.index)
        self.ModalBasis = pyfits.getdata(self.CDMS_BaseDir+str(self.CIAO_ID)+self.CDMS_ConfigDir+ 
                          'RecnOptimiser.ModalBasis.fits', ignore_missing_end=True)
        self.HOIM = pyfits.getdata(self.datadir+'/RecnOptimiser.HO_IM_%04d.fits'%self.index,
                                   ignore_missing_end=True)
        self.TT2HO = pyfits.getdata(self.CDMS_BaseDir+str(self.CIAO_ID)+self.CDMS_ConfigDir+
                                    'RecnOptimiser.TT2HO.fits', ignore_missing_end=True)
        self.iTT2HO = pyfits.getdata(self.CDMS_BaseDir+str(self.CIAO_ID)+self.CDMS_ConfigDir+
                                    'RecnOptimiser.ITT2HO.fits', ignore_missing_end=True)
        self.CM = pyfits.getdata(self.datadir+'/Recn.REC1.CM_%04d.fits'%self.index, 
                                 ignore_missing_end=True)
        self.CM /= self.controlGain
        self.CM[:60,:] += self.TT2HO.dot(self.CM[60:,:])
        self.TTM2Z = pyfits.getdata(self.CDMS_BaseDir+self.CDMS_ConfigDir+
                                    'RecnOptimiser.TTM2Z.fits', ignore_missing_end = True)
        self.Z2DM = pyfits.getdata(self.CDMS_BaseDir+str(self.CIAO_ID)+self.CDMS_ConfigDir+
                                   'RecnOptimiser.Z2DM.fits', ignore_missing_end=True)
        self.NZernikes = self.Z2DM.shape[1]
        self.ZIn = numpy.array([(i >= 3) & (i < 99) for i in range(self.NZernikes)])
        #self.DM2Z = linalg.pinv(self.Z2DM)
        self.DM2Z = pyfits.getdata(self.CDMS_BaseDir+str(self.CIAO_ID)+self.CDMS_ConfigDir+
                                   'RecnOptimiser.DM2Z.fits', ignore_missing_end=True)
        self.S2Z = self.DM2Z.dot(self.CM[:60])
        self.Z2S = linalg.pinv(self.S2Z)
        self.Voltage2Zernike = numpy.concatenate((self.DM2Z.T, self.TTM2Z.T))

    def extractModulation(self):
        self.modulatedFrames = numpy.append(numpy.diff(self.HODM, axis=0)[:,0] != 0, False)
        mean = numpy.mean(self.HODM[self.modulatedFrames==False], axis=0)
        self.modulatedHODM = self.HODM - mean
        mean = numpy.mean(self.Gradients[self.modulatedFrames==False], axis=0)
        self.modulatedGradients = self.Gradients - mean

    def extractIMsFromNoiseModulation(self, ax=None):
        PERIOD = self.loopTime[self.modulatedFrames==True][-1]/(len(self.modulatedFrames==True)-1)
        LOOPRATE = 1.0/PERIOD
        fourierHODM = numpy.fft.fftshift(numpy.fft.fft(self.modulatedHODM[self.modulatedFrames==True]))
        fourierGradients = numpy.fft.fftshift(numpy.fft.fft(self.modulatedGradients[self.modulatedFrames == True]))
        #fourierHODM = numpy.sqrt(fourierHODM * numpy.conj(fourierHODM))
        freq = numpy.fft.fftshift(numpy.fft.fftfreq(fourierHODM.shape[0], d=PERIOD))
        if ax != None:
            ax.clear()
            for H in fourierHODM.T:
                ax.plot(freq, H)
            ax.figure.show()
            raw_input()
            ax.clear()
            for G in fourierGradients.T:
                ax.plot(freq, G)
            ax.figure.show()
            raw_input()


        
class Anamorphose( object ):
    def __init__(self, directory='', RTC_Delay=0.5e-3, sqlCursor=None):
        self.dir = directory
        self.datadir = os.environ.get('CIAO_DATA')
        self.WFS_Frame = WFS_Frame()
        self.sqlCursor = sqlCursor

    def loadData(self, ax=None):
        self.header = pyfits.getheader(self.datadir+self.dir+'/CIAO_LOOP_0001.fits')
        self.CIAO_ID = self.header.get("ESO OCS SYS ID")
        self.nMeasurements = self.header.get('ESO TPL NEXP')
        self.Anamorphs = []
        self.IMs = []
        for i in range(self.nMeasurements):
            try:
                self.Anamorphs.append(Anamorph(directory=self.dir, index=i+1))
                self.Anamorphs[-1].loadData()
                self.Anamorphs[-1].extractModulation()
                self.IMs.append(self.Anamorphs[-1].extractIMsFromNoiseModulation(ax=ax))
            except:
                pass
        try:
            self.startTime = aptime.Time(self.header.get('ESO TPL START')).mjd
        except:
            headerTime = self.dir.split('/')
            HMS = headerTime[-1].split('-')[1]
            H = HMS[0:2]
            M = HMS[2:4]
            S = HMS[4:]
            self.startTime = aptime.Time(headerTime[2]+'T'+H+':'+M+':'+S).mjd
        """
        self.startTime = time.mktime(time.strptime(self.header.get('ESO TPL START'), 
                                     '%Y-%m-%dT%H:%M:%S'))
        self.time = loopData.field('Seconds')+data.field('USeconds')/100000.0
        self.time -= self.time[0]
        self.loopRate = self.header.get('ESO AOS LOOP RATE')
        self.controlGain = self.header.get('ESO AOS GLOBAL GAIN')
        #"""

    def addToDatabase(self):
        if self.sqlCursor == None:
            return
        sqlCommand = "SELECT * FROM CIAO_%d_Anamorphose WHERE TIMESTAMP = %.1f" % (self.CIAO_ID, 
                      self.startTime)
        self.sqlCursor.execute(sqlCommand)
        result = self.sqlCursor.fetchall()
        if len(result) == 0:
            values = {}
            values["TIMESTAMP"] = self.startTime
            values["DIRECTORY"] = self.dir
            values["CM_MODES"] = self.header.get('ESO AOS CM MODES CONTROLLED')
            values["WFS_GEOM"] = self.header.get('ESO AOS GEOM NAME')
            values["GAIN"] = self.header.get('ESO AOS GLOBAL GAIN')
            values["LOOPRATE"] = self.header.get('ESO AOS LOOP RATE')
            values["TTX_REFPOS"] = self.header.get('ESO AOS TTM REFPOS X')
            values["TTY_REFPOS"] = self.header.get('ESO AOS TTM REFPOS Y')
            values["WFS_MODE"] = self.header.get('ESO AOS WFS MODE')
            values["ALT"] = self.header.get('ESO TEL ALT')
            values["AZ"] = self.header.get('ESO TEL AZ')
            values["DROT_ENC"] = self.header.get('ESO INS DROT ENC')
            values["FILT_ENC"] = self.header.get('ESO INS FILT1 ENC')
            values["MSEL_ENC"] = self.header.get('ESO INS MSEL ENC')
            values["PMTIL_ENC"] = self.header.get('ESO INS PMTIL ENC')
            values["PMTIP_ENC"] = self.header.get('ESO INS PMTIP ENC')
            values["FLDLX"] = self.header.get('ESO INS FLDL X')
            values["FLDLY"] = self.header.get('ESO INS FLDL Y')
            values["FSM_A_U"] = self.header.get('ESO STS FSM1 GUIDE U')
            values["FSM_A_W"] = self.header.get('ESO STS FSM1 GUIDE W')
            values["FSM_A_X"] = self.header.get('ESO STS FSM1 POSX')
            values["FSM_A_Y"] = self.header.get('ESO STS FSM1 POSY')
            values["FSM_B_U"] = self.header.get('ESO STS FSM2 GUIDE U')
            values["FSM_B_W"] = self.header.get('ESO STS FSM2 GUIDE W')
            values["FSM_B_X"] = self.header.get('ESO STS FSM2 POSX')
            values["FSM_B_Y"] = self.header.get('ESO STS FSM2 POSY')
            values["VCM_A_U"] = self.header.get('ESO STS VCM1 GUIDE U')
            values["VCM_A_W"] = self.header.get('ESO STS VCM1 GUIDE W')
            values["VCM_A_X"] = self.header.get('ESO STS VCM1 POSX')
            values["VCM_A_Y"] = self.header.get('ESO STS VCM1 POSY')
            values["VCM_B_U"] = self.header.get('ESO STS VCM2 GUIDE U')
            values["VCM_B_W"] = self.header.get('ESO STS VCM2 GUIDE W')
            values["VCM_B_X"] = self.header.get('ESO STS VCM2 POSX')
            values["VCM_B_Y"] = self.header.get('ESO STS VCM2 POSY')
            values["M10_POSANG"] = self.header.get('ESO STS M10 POSANG')

            format_str= """INSERT INTO CIAO_"""+str(self.CIAO_ID)+"""_Anamorphose (TIMESTAMP, DIRECTORY, 
                 CM_MODES, WFS_GEOM, GAIN, LOOPRATE, TTX_REFPOS, TTY_REFPOS, WFS_MODE,
                 ALT, AZ, DROT_ENC, FILT_ENC, MSEL_ENC, PMTIL_ENC, PMTIP_ENC, FLDLX, FLDLY, FSM_A_U, 
                 FSM_A_W, FSM_A_X, FSM_A_Y, FSM_B_U, FSM_B_W, FSM_B_X, FSM_B_Y, VCM_A_U, VCM_A_W, 
                 VCM_A_X, VCM_A_Y, VCM_B_U, VCM_B_W, VCM_B_X, VCM_B_Y, M10_POSANG)
                 VALUES ('{timestamp}', '{directory}', '{cm_modes}', 
                 '{wfs_geom}', '{gain}', '{looprate}', '{ttx_refpos}', '{tty_refpos}', '{wfs_mode}', 
                 '{alt}', '{az}', '{derot_enc}', '{filt_enc}', '{msel_enc}', '{pmtil_enc}', 
                 '{pmtip_enc}', '{fldlx}', '{fldly}', '{fsm_a_u}', '{fsm_a_w}', '{fsm_a_x}','{fsm_a_y}',
                 '{fsm_b_u}', '{fsm_b_w}', '{fsm_b_x}', '{fsm_b_y}', '{vcm_a_u}',
                 '{vcm_a_w}', '{vcm_a_x}', '{vcm_a_y}', '{vcm_b_u}', '{vcm_b_w}',
                 '{vcm_b_x}', '{vcm_b_y}', '{m10}');"""
                
            sqlCommand = format_str.format(timestamp=values['TIMESTAMP'],directory=values['DIRECTORY'],
                 cm_modes=values['CM_MODES'], wfs_geom=values['WFS_GEOM'], gain=values['GAIN'],
                 looprate=values['LOOPRATE'], ttx_refpos=values['TTX_REFPOS'],
                 tty_refpos=values['TTY_REFPOS'], wfs_mode=values['WFS_MODE'],
                 alt=values['ALT'], az=values['AZ'], derot_enc=values['DROT_ENC'], 
                 filt_enc=values['FILT_ENC'], msel_enc=values['MSEL_ENC'],
                 pmtil_enc=values['PMTIL_ENC'], pmtip_enc=values['PMTIP_ENC'],
                 fldlx=values['FLDLX'], fldly=values['FLDLY'], fsm_a_u=values['FSM_A_U'],
                 fsm_a_w=values['FSM_A_W'], fsm_a_x=values['FSM_A_X'],
                 fsm_a_y=values['FSM_A_Y'], fsm_b_u=values['FSM_B_U'],
                 fsm_b_w=values['FSM_B_W'], fsm_b_x=values['FSM_B_X'],
                 fsm_b_y=values['FSM_B_Y'], vcm_a_u=values['VCM_A_U'],
                 vcm_a_w=values['VCM_A_W'], vcm_a_x=values['VCM_A_X'],
                 vcm_a_y=values['VCM_A_Y'], vcm_b_u=values['VCM_B_U'],
                 vcm_b_w=values['VCM_B_W'], vcm_b_x=values['VCM_B_X'],
                 vcm_b_y=values['VCM_B_Y'], m10=values['M10_POSANG'])
            
            self.sqlCursor.execute(sqlCommand)
 

class AVC( object ):
    def __init__(self, PafFile='', AVCFile=''):
        self.PafFile = PafFile
        self.AVCFile = AVCFile
        if PafFile != '':
            self.parsePaf()
        if AVCFile != '':
            self.parseAVC()
        

    def parsePaf(self):
        f = open(self.PafFile, 'r')
        for line in f.readlines():
            if 'HOCtrAVC.CFG.DYNAMIC.MODE' in line:
                print 'Hi!'

    def parseAVC(self):
        data = pyfits.getdata(self.AVCFile)
        self.FC = data.field('FrameCounter')
        self.X1 = data.field('X1')
        self.X2 = data.field('X2')
        self.X3 = data.field('X3')
        self.X4 = data.field('X4')
        self.Omega = data.field('Omega')
        self.freq = self.Omega/(2.0*3.14159)
        self.FilterInput = data.field('FilterInput')
        self.FilterOutput = data.field('FilterOutput')

        self.Modes = {}
        for i in range(5):
            self.Modes[i] = {}
            self.Modes[i]['Frequencies'] = self.freq[0,i*10:(i+1)*10]
            enabled = []
            tripped = []
            for j in range(10):
                power = numpy.unique(self.FilterOutput[:,i*10+j])
                enabled.append(len(power) > 1)
                tripped.append((len(power) > 1) and (0 in power))
            self.Modes[i]['Enabled'] = numpy.array(enabled)
            self.Modes[i]['Tripped'] = numpy.array(tripped)

    def setModes(self, modes):
        self.Modes = {}
        for i in range(len(modes)):
            self.Modes[i] = {}
            self.Modes[i]['Frequencies'] = range(10)
            self.Modes[i]['Enabled'] = numpy.array([True for j in range(10)])
            self.Modes[i]['Tripped'] = numpy.array([False for j in range(10)])

class DataLogger( object ):
    def __init__(self, directory='', CDMS_BaseDir='', CDMS_ConfigDir='', RTC_Delay=0.5e-3, 
                 sqlCursor=None):
        self.dir = directory
        self.datadir = os.environ.get('CIAO_DATA')
        self.CDMS_BaseDir = os.environ.get('CDMS_BASEDIR')
        self.CDMS_ConfigDir = os.environ.get('CDMS_CONFIGDIR')
        self.MATLAB_CIAO_DATABASE = os.environ.get('MATLAB_CIAO_DATABASE')
        #self.CDMS_BaseDir = CDMS_BaseDir
        #self.CDMS_ConfigDir = CDMS_ConfigDir
        self.RTC_Delay = RTC_Delay
        self.ApertureDiameter = 8.0
        self.WFS_Frame = WFS_Frame()
        self.sqlCursor = sqlCursor

    def loadData(self):
        self.header = pyfits.getheader(self.datadir+self.dir+'/CIAO_LOOP_0001.fits')
        self.CIAO_ID = self.header.get("ESO OCS SYS ID")
        self.cimdatDir = self.MATLAB_CIAO_DATABASE+'UT%d/OFFAXIS/Matrices/Laptop/' % self.CIAO_ID
        self.LambdaSeeing = 0.5e-6
        self.LambdaStrehl = 2.2e-6
        self.Arcsec = 180.0*3600.0/numpy.pi
        self.ReferenceSeeing = 1.0/self.Arcsec
        self.L0 = numpy.inf
        self.ApertureDiameter = 8.0
        self.r0 = self.LambdaSeeing/self.ReferenceSeeing
        data = pyfits.getdata(self.datadir+self.dir+'/CIAO_LOOP_0001.fits')
        self.Intensities = data.field('Intensities')
        self.Gradients = data.field('Gradients')
        self.HODM = data.field('HODM_Positions')
        self.TTM = data.field('ITTM_Positions')
        self.FrameCounter = data.field('FrameCounter')
        self.AVC = AVC(AVCFile=self.datadir+self.dir+'/CIAO_AVC_0001.fits')
        try:
            self.startTime = aptime.Time(self.header.get('ESO TPL START')).mjd
        except:
            headerTime = self.dir.split('/')
            HMS = headerTime[-1].split('-')[1]
            H = HMS[0:2]
            M = HMS[2:4]
            S = HMS[4:]
            self.startTime = aptime.Time(headerTime[2]+'T'+H+':'+M+':'+S).mjd
        self.time = data.field('Seconds')+data.field('USeconds')/100000.0
        self.time -= self.time[0]
        self.loopRate = self.header.get('ESO AOS LOOP RATE')
        self.controlGain = self.header.get('ESO AOS GLOBAL GAIN')
        try:
            self.S2M = pyfits.getdata(self.datadir+self.dir+'/RecnOptimiser.S2M_0001.fits')
        except:
            self.S2M = pyfits.getdata(self.datadir+self.dir+'/CLMatrixOptimiser.S2M_0001.fits')
        self.ModalBasis = pyfits.getdata(self.CDMS_BaseDir+str(self.CIAO_ID)+self.CDMS_ConfigDir+ 
                          'RecnOptimiser.ModalBasis.fits', ignore_missing_end=True)
        try:
            self.HOIM = pyfits.getdata(self.datadir+self.dir+'/RecnOptimiser.HO_IM_0001.fits',
                                   ignore_missing_end=True)
        except:
            self.HOIM = pyfits.getdata(self.datadir+self.dir+'/CLMatrixOptimiser.HO_IM_0001.fits',
                                   ignore_missing_end=True)
        self.TT2HO = pyfits.getdata(self.CDMS_BaseDir+str(self.CIAO_ID)+self.CDMS_ConfigDir+
                                    'RecnOptimiser.TT2HO.fits', ignore_missing_end=True)
        self.iTT2HO = pyfits.getdata(self.CDMS_BaseDir+str(self.CIAO_ID)+self.CDMS_ConfigDir+
                                    'RecnOptimiser.ITT2HO.fits', ignore_missing_end=True)
        self.CM = pyfits.getdata(self.datadir+self.dir+'/Recn.REC1.CM_0001.fits', 
                                 ignore_missing_end=True)
        if self.controlGain > 0.0:
            self.CM /= self.controlGain
        self.CM[:60,:] += self.TT2HO.dot(self.CM[60:,:])
        self.TTM2Z = pyfits.getdata(self.cimdatDir+
                                    'cimdatTTM2Z.fits', ignore_missing_end = True)
        self.Z2DM = pyfits.getdata(self.CDMS_BaseDir+str(self.CIAO_ID)+self.CDMS_ConfigDir+
                                   'RecnOptimiser.Z2DM.fits',
                                   ignore_missing_end=True)[:,:60]
        self.NZernikes = self.Z2DM.shape[1]
        self.ZIn = numpy.array([(i >= 3) & (i < 99) for i in range(self.NZernikes)])
        #self.ZIn = numpy.array([(i >= 3) & (i < 99) for i in range(59)])
        #self.DM2Z = linalg.pinv(self.Z2DM)
        self.DM2Z = pyfits.getdata(self.CDMS_BaseDir+str(self.CIAO_ID)+self.CDMS_ConfigDir+
                                   'RecnOptimiser.DM2Z.fits', ignore_missing_end=True)
        #self.S2Z = self.DM2Z.dot(self.CM[:60])
        #self.Z2S = linalg.pinv(self.S2Z)
        self.Voltage2Zernike = numpy.concatenate((self.DM2Z.T,
            self.TTM2Z.T))[:,:59]
        self.Zernike2Slopes =  pyfits.getdata(self.cimdatDir+'cimdatZernike2Slopes.fits')
        self.S2Z = linalg.pinv(self.Zernike2Slopes)
        try:
            self.REFSLP = pyfits.getdata(self.datadir+self.dir+'/Acq.DET1.REFSLP_0001.fits')
        except:
            self.REFSLP = None

    def __lt__(self, other):
        if hasattr(other, startTime):
            return self.startTime < other.startTime

    def __le__(self, other):
        if hasattr(other, startTime):
            return self.startTime <= other.startTime

    def __gt__(self, other):
        if hasattr(other, startTime):
            return self.startTime > other.startTime

    def __ge__(self, other):
        if hasattr(other, startTime):
            return self.startTime >= other.startTime

    def __eq__(self, other):
        if hasattr(other, startTime):
            return self.startTime == other.startTime

    def __ne__(self, other):
        if hasattr(other, startTime):
            return self.startTime != other.startTime

    def getRefSlopeZernikes(self):
        #self.rotatedZernikes = self.S2Z.dot(self.REFSLP.T)
        self.rotatedZernikes = self.REFSLP.dot(self.S2Z)
        
        return #self.S2Z.dot(self.REFSLP.T)

    def updateTimeStamp(self, oldTimeStamp):
        sqlCommand = "UPDATE CIAO_%d_DataLoggers SET TIMESTAMP = %.7f WHERE TIMESTAMP =%.1f" % (self.CIAO_ID, self.startTime, oldTimeStamp)
        self.sqlCursor.execute(sqlCommand)

    def addToDatabase(self):
        if self.sqlCursor == None:
            return
        sqlCommand = "SELECT * FROM CIAO_%d_DataLoggers WHERE TIMESTAMP = %.1f" % (self.CIAO_ID, 
                      self.startTime)
        self.sqlCursor.execute(sqlCommand)
        result = self.sqlCursor.fetchall()
        if len(result) == 0:
            values = {}
            values['timestamp'] = self.startTime
            values['directory'] = self.dir
            values['asm_seeing'] = self.header.get('ESO TEL AMBI FWHM')
            values['avc_state'] = self.header.get('ESO AOS AVC LOOP ST') == True
            values['cm_modes'] = self.header.get('ESO AOS CM MODES CONTROLLED')
            values['wfs_geom'] = self.header.get('ESO AOS GEOM NAME')
            values['gain'] = self.header.get('ESO AOS GLOBAL GAIN')
            values['awf_enable'] = self.header.get('ESO AOS HOCTR AWF ENABLE')
            values['garbage_gain'] = self.header.get('ESO AOS HOCTR GARBAGE GAIN')
            values['ki'] = self.header.get('ESO AOS HOCTR KI')
            values['kt'] = self.header.get('ESO AOS HOCTR KT')
            values['pra_enable'] = self.header.get('ESO AOS HOCTR PRA ENABLE')
            values['pra_gain'] = self.header.get('ESO AOS HOCTR PRA GAIN')
            values['sma_enable'] = self.header.get('ESO AOS HOCTR SMA ENABLE')
            values['sma_high'] = self.header.get('ESO AOS HOCTR SMA HIGH')
            values['sma_iterations'] = self.header.get('ESO AOS HOCTR SMA ITERATIONS')
            values['sma_low'] = self.header.get('ESO AOS HOCTR SMA LOW')
            values['tt_ki'] = self.header.get('ESO AOS HOCTR TT KI')
            values['tt_kt'] = self.header.get('ESO AOS HOCTR TT KT')
            values['looprate'] = self.header.get('ESO AOS LOOP RATE')
            values['tt_loop_state'] = self.header.get('ESO AOS TT LOOP ST')
            values['ttx_refpos'] = self.header.get('ESO AOS TTM REFPOS X')
            values['tty_refpos'] = self.header.get('ESO AOS TTM REFPOS Y')
            values['vib_sr'] = self.header.get('ESO AOS VIB SR')
            values['wfs_mode'] = self.header.get('ESO AOS WFS MODE')
            values['im_trk_mode'] = self.header.get('ESO OCS IM TRK MODE')
            values['alt'] = self.header.get('ESO TEL ALT')
            values['az'] = self.header.get('ESO TEL AZ')
            values['dit'] = self.header.get('ESO DET DIT')
            values['derot_enc'] = self.header.get('ESO INS DROT ENC')
            values['filt_enc'] = self.header.get('ESO INS FILT1 ENC')
            values['msel_enc'] = self.header.get('ESO INS MSEL ENC')
            values['msel_name'] = self.header.get('ESO INS MSEL NAME')
            values['pmtil_enc'] = self.header.get('ESO INS PMTIL ENC')
            values['pmtip_enc'] = self.header.get('ESO INS PMTIP ENC')
            values['fldlx'] = self.header.get('ESO INS FLDL X')
            values['fldly'] = self.header.get('ESO INS FLDL Y')
            values['fsm_a_u'] = self.header.get('ESO STS FSM1 GUIDE U')
            values['fsm_a_w'] = self.header.get('ESO STS FSM1 GUIDE W')
            values['fsm_a_x'] = self.header.get('ESO STS FSM1 POSX')
            values['fsm_a_y'] = self.header.get('ESO STS FSM1 POSY')
            values['fsm_b_u'] = self.header.get('ESO STS FSM2 GUIDE U')
            values['fsm_b_w'] = self.header.get('ESO STS FSM2 GUIDE W')
            values['fsm_b_x'] = self.header.get('ESO STS FSM2 POSX')
            values['fsm_b_y'] = self.header.get('ESO STS FSM2 POSY')
            values['vcm_a_u'] = self.header.get('ESO STS VCM1 GUIDE U')
            values['vcm_a_w'] = self.header.get('ESO STS VCM1 GUIDE W')
            values['vcm_a_x'] = self.header.get('ESO STS VCM1 POSX')
            values['vcm_a_y'] = self.header.get('ESO STS VCM1 POSY')
            values['vcm_b_u'] = self.header.get('ESO STS VCM2 GUIDE U')
            values['vcm_b_w'] = self.header.get('ESO STS VCM2 GUIDE W')
            values['vcm_b_x'] = self.header.get('ESO STS VCM2 POSX')
            values['vcm_b_y'] = self.header.get('ESO STS VCM2 POSY')
            values['m10'] = self.header.get('ESO STS M10 POSANG')
            values['rhum'] = self.header.get('ESO TEL AMBI RHUM')
            values['temp'] = self.header.get('ESO TEL AMBI TEMP')
            values['theta0'] = self.header.get('ESO TEL AMBI THETA0')
            values['winddir'] = self.header.get('ESO TEL AMBI WINDDIR')
            values['windsp'] = self.header.get('ESO TEL AMBI WINDSP')
            values['prltic'] = self.header.get('ESO TEL PRLTIC')
            values['pup_trk'] = self.header.get('ESO PUP TRK ST')
                
            format_str= """INSERT INTO CIAO_"""+str(self.CIAO_ID)+"""_DataLoggers (TIMESTAMP, DIRECTORY, 
                 ASM_SEEING, AVC_STATE, CM_MODES, WFS_GEOM, GAIN,
                 HOCTR_AWF_ENABLE, HOCTR_GARBAGE_GAIN, HOCTR_KI, HOCTR_KT, 
                 HOCTR_PRA_ENABLE, HOCTR_PRA_GAIN, HOCTR_SMA_ENABLE, HOCTR_SMA_HIGH,
                 HOCTR_SMA_ITERATIONS, HOCTR_SMA_LOW, HOCTR_TT_KI, HOCTR_TT_KT,
                 LOOPRATE, TT_LOOP_STATE, TTX_REFPOS, TTY_REFPOS, VIB_SR, WFS_MODE,
                 IM_TRK_MODE, ALT, AZ, DIT, DROT_ENC, FILT_ENC, MSEL_ENC, MSEL_NAME,
                 PMTIL_ENC, PMTIP_ENC, FLDLX, FLDLY, FSM_A_U, FSM_A_W, FSM_A_X, FSM_A_Y,
                 FSM_B_U, FSM_B_W, FSM_B_X, FSM_B_Y, VCM_A_U, VCM_A_W, VCM_A_X, VCM_A_Y,
                 VCM_B_U, VCM_B_W, VCM_B_X, VCM_B_Y, M10_POSANG, RHUM, TEMP, THETA0, 
                 WINDDIR, WINDSP, PRLTIC, PUP_TRK)
                 VALUES ('{timestamp}', '{directory}', '{asm_seeing}',
                 '{avc_state}', '{cm_modes}', 
                 '{wfs_geom}', '{gain}', '{awf_enable}', '{garbage_gain}', '{ki}',
                 '{kt}', '{pra_enable}', '{pra_gain}', '{sma_enable}', '{sma_high}',
                 '{sma_iterations}', '{sma_low}', '{tt_ki}', '{tt_kt}', '{looprate}',
                 '{tt_loop_state}', '{ttx_refpos}', '{tty_refpos}', '{vib_sr}',
                 '{wfs_mode}', '{im_trk_mode}', '{alt}', '{az}', '{dit}', '{derot_enc}',
                 '{filt_enc}', '{msel_enc}', '{msel_name}', '{pmtil_enc}', '{pmtip_enc}',
                 '{fldlx}', '{fldly}', '{fsm_a_u}', '{fsm_a_w}', '{fsm_a_x}', '{fsm_a_y}',
                 '{fsm_b_u}', '{fsm_b_w}', '{fsm_b_x}', '{fsm_b_y}', '{vcm_a_u}',
                 '{vcm_a_w}', '{vcm_a_x}', '{vcm_a_y}', '{vcm_b_u}', '{vcm_b_w}',
                 '{vcm_b_x}', '{vcm_b_y}', '{m10}', '{rhum}', '{temp}', '{theta0}',
                 '{winddir}', '{windsp}', '{prltic}', '{pup_trk}');"""
                
            sqlCommand = format_str.format(timestamp=values['timestamp'],directory=values['directory'],
                 asm_seeing=values['asm_seeing'], avc_state=values['avc_state'], 
                 cm_modes=values['cm_modes'],
                 wfs_geom=values['wfs_geom'], gain=values['gain'],
                 awf_enable=values['awf_enable'], garbage_gain=values['garbage_gain'],
                 ki=values['ki'], kt=values['kt'], pra_enable=values['pra_enable'],
                 pra_gain=values['pra_gain'], sma_enable=values['sma_enable'],
                 sma_high=values['sma_high'], sma_iterations=values['sma_iterations'],
                 sma_low=values['sma_low'], tt_ki=values['tt_ki'],
                 tt_kt=values['tt_kt'], looprate=values['looprate'],
                 tt_loop_state=values['tt_loop_state'], ttx_refpos=values['ttx_refpos'],
                 tty_refpos=values['tty_refpos'], vib_sr=values['vib_sr'],
                 wfs_mode=values['wfs_mode'], im_trk_mode=values['im_trk_mode'],
                 alt=values['alt'], az=values['az'], dit=values['dit'],
                 derot_enc=values['derot_enc'], filt_enc=values['filt_enc'],
                 msel_enc=values['msel_enc'], msel_name=values['msel_name'],
                 pmtil_enc=values['pmtil_enc'], pmtip_enc=values['pmtip_enc'],
                 fldlx=values['fldlx'], fldly=values['fldly'], fsm_a_u=values['fsm_a_u'],
                 fsm_a_w=values['fsm_a_w'], fsm_a_x=values['fsm_a_x'],
                 fsm_a_y=values['fsm_a_y'], fsm_b_u=values['fsm_b_u'],
                 fsm_b_w=values['fsm_b_w'], fsm_b_x=values['fsm_b_x'],
                 fsm_b_y=values['fsm_b_y'], vcm_a_u=values['vcm_a_u'],
                 vcm_a_w=values['vcm_a_w'], vcm_a_x=values['vcm_a_x'],
                 vcm_a_y=values['vcm_a_y'], vcm_b_u=values['vcm_b_u'],
                 vcm_b_w=values['vcm_b_w'], vcm_b_x=values['vcm_b_x'],
                 vcm_b_y=values['vcm_b_y'], m10=values['m10'], rhum=values['rhum'],
                 temp=values['temp'], theta0=values['theta0'], winddir=values['winddir'],
                 windsp=values['windsp'], prltic=values['prltic'], pup_trk=values['pup_trk'])
            
            self.sqlCursor.execute(sqlCommand)

    def calculatePhotometricPupilOffset(self):
        self.WFS_Frame.generateIntensityImage(intensities = numpy.mean(self.Intensities, axis=0))
        self.innerRing = numpy.array(center_of_mass(self.WFS_Frame.innerRingImage)) - \
                              numpy.array([4.0, 4.0])
        self.outerRing = numpy.array(center_of_mass(self.WFS_Frame.outerRingImage)) - \
                              numpy.array([4.0, 4.0])

    def computePSFPerformance(self):
        self.FrequencyResolution = 0.5
        self.synchronizeData()
        self.zernikeSpace()
        #self.Piston = self.self.HODM    - don't have piston matrices loaded in yet...


    def computeStrehlSegments(self, DITTimes=[]):
        t1 = self.startTime
        t2 = self.startTime
        oldGradients = self.Gradients
        oldTTM = self.TTM
        oldHODM = self.HODM
        strehl = []
        seeing = []
        for t in DITTimes:
            t1 = t2
            t2 = t
            self.TTM = oldTTM
            self.HODM = oldHODM
            self.Gradients = oldGradients
            self.clipData(startTime=t1, stopTime=t2)
            if len(self.TTM) > 2000:
                print len(self.TTM)
                self.synchronizeData()
                self.zernikeSpace()
                self.computePSD(source='ZSlopes')
                self.computePSD(source='ZCommands')
                self.AORejectionFunction()
                self.combinePSDs()
                self.computeKolmogorovCovar()
                self.zernikePropagation()
                self.noiseEvaluation()
                self.seeingEstimation()
                self.computeSpectralSlope((3, 15))
                self.estimateStrehlRatio()
                strehl.append(self.Strehl)
                seeing.append(self.Seeing*self.Arcsec)
        return numpy.array(strehl), numpy.array(seeing)

    def clipData(self, startTime=0.0, stopTime=0.0):
        frames = (self.time+self.startTime > startTime) & (self.time+self.startTime < stopTime)
        self.TTM = self.TTM[frames, :]
        self.Gradients = self.Gradients[frames,:]
        self.HODM = self.HODM[frames, :]

    def computeStrehl(self, saveData=False, skipLongRecords=False):
        if (skipLongRecords and(self.HODM.shape[0] > 30000)):
            self.Seeing = numpy.nan
            self.Strehl = numpy.nan
            self.Arcsec = 1.0
            print "Skipping long record!"
            return
        self.synchronizeData()
        self.zernikeSpace()
        self.computePSD(source='ZSlopes')
        self.computePSD(source='ZCommands')
        self.AORejectionFunction()
        self.combinePSDs()
        self.computeKolmogorovCovar()
        self.zernikePropagation()
        self.noiseEvaluation()
        self.seeingEstimation()
        self.computeSpectralSlope((3,15))
        self.estimateStrehlRatio()

        if (saveData and numpy.isfinite(self.Strehl) and numpy.isfinite(self.Seeing) and
                     numpy.isfinite(self.Tau0)):
            sqlCommand = "UPDATE CIAO_%d_DataLoggers SET STREHL = %.4f, SEEING = %.4f, TAU0 = %.4f, TERR = %.4f WHERE TIMESTAMP = %.7f;" % (self.CIAO_ID, self.Strehl, self.Seeing*self.Arcsec, self.Tau0, self.TemporalError, self.startTime)
            self.sqlCursor.execute(sqlCommand)
    
    def measureVibs(self, frequencies=[], modes='AVC', saveData=False):
        self.vibPower = {}
        if modes=='AVC':
            extractedModes = [0, 1, 2, 7, 8]
        elif modes=='ALL':
            extractedModes = range(self.ZPowerCommands.shape[0])
            self.AVC.setModes(extractedModes)
        elif modes=='TEST':
            extractedModes = [1, 2]
            self.AVC.setModes(extractedModes)
        for Mode, i in zip(self.AVC.Modes.itervalues(), extractedModes):
            self.vibPower[i] = {}
            commPower = {}
            slopesPower = {}
            #self.flattenPSD(i)
            #for j in range(10):
            #    if Mode['Enabled'][j]:
            #        freq = Mode['Frequencies'][j]
            #        commPower[freq]=self.computeCommVibPower(i, Mode['Frequencies'][j])
            #        slopesPower[freq]=self.computeSlopesVibPower(i, Mode['Frequencies'][j])
            for f in frequencies:
                commPower[f] = self.computeCommVibPower(i, f)
                slopesPower[f] = self.computeSlopesVibPower(i, f)
            self.vibPower[i]['CommPower'] = commPower
            self.vibPower[i]['SlopesPower'] = slopesPower

    def flattenPSD(self, i):
        window = self.ZPowerFrequencies > 15.0

        SlopesA = numpy.log10(self.ZPowerSlopes[i,window])
        CommsA = numpy.log10(self.ZPowerCommands[i,window])
        model = numpy.array([numpy.ones(self.ZPowerFrequencies[window].shape), numpy.log10(self.ZPowerFrequencies[window])])
        im = numpy.linalg.pinv(model)
        SlopeSpectralSlope = im.T.dot(SlopesA.T)
        CommSpectralSlope = im.T.dot(CommsA.T)

        fullModel = numpy.array([numpy.ones(self.ZPowerFrequencies.shape), numpy.log10(self.ZPowerFrequencies)])
        self.ZPowerSlopes[i] = self.ZPowerSlopes[i]/10.0**SlopeSpectralSlope.dot(fullModel)
        self.ZPowerCommands[i] = self.ZPowerCommands[i]/10.0**CommSpectralSlope.dot(fullModel)


    def computeCommVibPower(self, mode, freq):
        Window = (self.ZPowerFrequencies > (freq - 1.0)) & (self.ZPowerFrequencies < 
                  (freq + 1.0))
        retval = scipy.integrate.trapz(self.ZPowerCommands[mode, Window],
                 x=self.ZPowerFrequencies[Window])
        return retval

    def computeSlopesVibPower(self, mode, freq):
        Window = (self.ZPowerFrequencies > (freq - 1.0)) & (self.ZPowerFrequencies < 
                  (freq + 1.0))
        retval = scipy.integrate.trapz(self.ZPowerSlopes[mode, Window],
                 x=self.ZPowerFrequencies[Window])
        return retval


    def computeTTPS(self, ax = None, freq = 2.2, saveData=False,
            returnTipTilt=False, Tip=None, Tilt=None):
        if Tip != None:
            A = numpy.log10(numpy.sum(numpy.array([10.0**Tip, 10.0**Tilt]), axis=0))
            m = numpy.array([numpy.ones(self.ZPowerFrequencies.shape),
                numpy.log10(self.ZPowerFrequencies)])
            im = numpy.linalg.pinv(m)
            spectralSlopes = im.T.dot(A.T)

        else:
            A = numpy.log10(numpy.sum(self.ZPowerCommands[0:2,], axis=0))
            m = numpy.array([numpy.ones(self.ZPowerFrequencies.shape),
                numpy.log10(self.ZPowerFrequencies)])
            im = numpy.linalg.pinv(m)
            spectralSlopes = im.T.dot(A.T)

        TipModel = spectralSlopes[0] + numpy.log10(self.ZPowerFrequencies)*spectralSlopes[1]
        TiltModel = spectralSlopes[0] + numpy.log10(self.ZPowerFrequencies)*spectralSlopes[1]
        if ax != None:
            ax.clear()
            ax.plot(numpy.log10(self.ZPowerFrequencies),
                    numpy.log10(self.ZPowerCommands[0]) - TipModel, color = 'r')
            ax.plot(numpy.log10(self.ZPowerFrequencies),
                    numpy.log10(self.ZPowerCommands[1]) - TiltModel, color = 'b')
            ax.plot(numpy.log10([freq, freq]), [-1, 1], color = 'g')
            ax.plot(numpy.log10([freq, freq]) - 0.15, [-1, 1], color = 'g')
            ax.plot(numpy.log10([freq, freq]) + 0.15, [-1, 1], color = 'g')
                
            ax.figure.show()
            raw_input()

        logFreq = numpy.log10(self.ZPowerFrequencies)
        Window = (logFreq > (numpy.log10(freq) - 0.15)) & (logFreq <
                (numpy.log10(freq) + 0.15))

        # Power in the Tip command
        logTip = numpy.log10(self.ZPowerCommands[0]) - TipModel
        logTip -= numpy.median(logTip)
        TipPower = scipy.integrate.trapz(logTip[Window], x=10.0**logFreq[Window])
        
        # Power in the Tilt command
        logTilt = numpy.log10(self.ZPowerCommands[1]) - TiltModel
        logTilt -= numpy.median(logTilt)
        TiltPower = scipy.integrate.trapz(logTilt[Window], x=10.0**logFreq[Window])

        print self.dir, TipPower, TiltPower
        if (saveData and numpy.isfinite(TipPower) and numpy.isfinite(TiltPower)):
            sqlCommand = "UPDATE CIAO_%d_DataLoggers SET TIPPOWER = %.2f, TILTPOWER = %.2f WHERE TIMESTAMP = %.7f;" % (self.CIAO_ID, TipPower, TiltPower, self.startTime)
            self.sqlCursor.execute(sqlCommand)

        if returnTipTilt:
            return numpy.log10(self.ZPowerCommands[0]), numpy.log10(self.ZPowerCommands[1]), numpy.log10(self.ZPowerFrequencies)


    def synchronizeData(self):
        self.TTM = self.TTM[:-4,:]
        self.HODM = self.HODM[2:-2,:]
        self.Gradients = self.Gradients[3:-1,:]
        
        self.TTM -= numpy.mean(self.TTM, axis=0)
        self.HODM -= numpy.mean(self.HODM, axis=0)
        self.Gradients -= numpy.mean(self.Gradients, axis=0)
        
    def zernikeSpace(self):
        self.ZSlopes = self.Gradients.dot(self.S2Z)
        self.ZCommands = numpy.concatenate((self.HODM.T, self.TTM.T)).T.dot(self.Voltage2Zernike)
        
    def pupilIllumination(self, ax):
        Top = [0, 1, 2, 3, 4]
        Bottom = [63, 64, 65, 66, 67]
        Left = [12, 21, 30, 39, 48]
        Right = [20, 29, 38, 47, 56]
        NFrames = self.Intensities.shape[0]

        for side in zip([Top, Bottom, Left, Right]):
            Illumination = numpy.average(self.Intensities[:, side], axis=2)
            Illumination = FourierDetrend(Illumination, 'TwoPoints')
            Power = numpy.abs(numpy.fft.fft(Illumination/NFrames))**2.0
            Freq = numpy.fft.fftfreq(NFrames, d=1.0/self.loopRate)
            Power = numpy.fft.fftshift(Power)
            Freq = numpy.fft.fftshift(Freq)
            Power =  2.0*Power[Freq >= 0]
            Freq = Freq[Freq>=0]
            Power /= 2.0
            Power = Power /numpy.mean(numpy.diff(Freq))
            ax.plot(Freq, Power)

        ax.set_yscale('log')
        #ax.figure.show()

    def computePSD(self, source):
        NWindows = 1
        Resolution = 0.5
        if source == 'ZSlopes':
            data = self.ZSlopes
        elif source == 'ZCommands':
            data = self.ZCommands
        NWindows = computeNWindows(data, self.loopRate, Resolution)
        NSeries = data.shape[1]
        data = FourierDetrend(data, 'TwoPoints') 

        NFrames = data.shape[0]
        Power = numpy.abs(numpy.fft.fft(data.T/NFrames))**2.0
        Freq = numpy.fft.fftfreq(NFrames, d=1.0/self.loopRate)

        Power = numpy.fft.fftshift(Power, axes=1)
        Freq = numpy.fft.fftshift(Freq)
        Power = 2.0*Power[:,Freq >= 0]
        Freq = Freq[Freq >= 0]
        Power[:,0] /= 2.0
        Power = Power / numpy.mean(numpy.diff(Freq))

        NSamplesPerWindow = numpy.floor(Power.shape[1]/NWindows)
        TotalSampleNumber = NWindows * NSamplesPerWindow
        Power = Power.T

        if NWindows > 1:
            Power = Power[:TotalSampleNumber, :]
            Power = Power.reshape([NWindows, NSamplesPerWindow, NSeries], order='F')
            Power = numpy.mean(Power, axis=0)

            Freq = Freq[:TotalSampleNumber]
            Freq = Freq.reshape([NWindows, NSamplesPerWindow], order = 'F')
            Freq = numpy.mean(Freq, axis = 0)

        if source == 'ZSlopes':
            self.ZPowerSlopes = Power
            self.ZPowerFrequencies = Freq
            self.ZPowerdFrequencies = numpy.mean(numpy.diff(self.ZPowerFrequencies))
        elif source == 'ZCommands':
            self.ZPowerCommands = Power
        return

    def AORejectionFunction(self):     # AORejection, WFS, RTC, NoisePropagation
        self.WFS2 = numpy.sum(self.ZPowerSlopes) * self.ZPowerdFrequencies
        Omega = 2*numpy.pi*self.ZPowerFrequencies
        iOmega = (0.0 + 1.0j)*Omega
        
        RTC_Period = 1.0/self.loopRate
        CycleDelay = numpy.exp(-RTC_Period*iOmega)
        CPU_Delay = numpy.exp(-self.RTC_Delay*iOmega)

        GainNumerator = numpy.array([self.controlGain])
        GainDenominator = numpy.array([1, -1])
        Numerator = GetTF(GainNumerator, CycleDelay)
        Denominator = GetTF(GainDenominator, CycleDelay)
        Controller = Numerator/Denominator


        DM = 1

        Ratio = self.RTC_Delay/RTC_Period
        N = numpy.floor(Ratio)
        alpha = Ratio - N

        RSinPhi = -alpha*numpy.sin(Omega*RTC_Period)
        RCosPhi = 1.0-alpha*(1.0-numpy.cos(Omega*RTC_Period))
        Phi = numpy.arctan2(RSinPhi, RCosPhi)
        R = numpy.sqrt(RSinPhi**2.0 + RCosPhi**2.0)

        self.WFS = R*numpy.exp((0.0+1.0j)*Phi)*CycleDelay**(N+1)

        self.RTC = Controller * DM
        OpenLoop_TF = self.WFS*self.RTC
        self.AORejection = 1.0/(1.0+OpenLoop_TF)
        self.NoisePropagation = self.RTC*self.AORejection

        self.AORejection[numpy.isnan(self.AORejection)]= 0.0
        self.WFS[numpy.isnan(self.WFS)] = 0.0
        self.RTC[numpy.isnan(self.RTC)] = 0.0
        self.NoisePropagation[numpy.isnan(self.NoisePropagation)] = 0.0
        self.AORejection = numpy.abs(self.AORejection)**2.0

    #def measureVibrations(self)

    def combinePSDs(self):
        RTC2 = numpy.abs(self.RTC)**2.0

        RTC2S = numpy.array([RTC2 * self.ZPowerSlopes[:,n] for n in range(self.ZPowerSlopes.shape[1])]).T
        self.ZPowerCommands = (RTC2S.T + [RTC2*self.ZPowerCommands[:,n] for n in
            range(self.ZPowerCommands.shape[1])])/(1.0+RTC2).T
        self.ZPowerSlopes = self.ZPowerCommands/RTC2
        

    def computeKolmogorovCovar(self):
        D_L0 = self.ApertureDiameter/self.L0

        infOuterScaleCovar, radial_orders = newKTaiaj(self.NZernikes+1, self.NZernikes+1, self.ApertureDiameter, self.r0)
        infOuterScaleCovar[:,0] = 0.0
        infOuterScaleCovar[0,:] = 0.0

        max_radial_order = numpy.int(numpy.max(radial_orders))

        reduc_var = OuterScale2VarianceReduction(D_L0, max_radial_order)
        reduction_factors = numpy.ones(radial_orders.shape)

        interp = scipy.interpolate.interp1d(numpy.array(range(max_radial_order))+1, reduc_var)
        blah = radial_orders != 0.0
        reduction_factors[blah] = interp(radial_orders[blah])

        self.KolmogorovCovar = infOuterScaleCovar * reduction_factors
        FullVariance, n = KolmogorovVariance(self.ApertureDiameter, self.L0, self.r0, 1000)

        self.FullVariance = numpy.sum(FullVariance)
        self.FractionOfTotalVariance = numpy.real(numpy.sum(numpy.diag(self.KolmogorovCovar))/self.FullVariance)

        self.OffControl = numpy.sum(numpy.diag(self.KolmogorovCovar))*(1.0/self.FractionOfTotalVariance -1.0)
        self.KolmogorovCovar = self.KolmogorovCovar[1:, 1:]

        self.KolmogorovCovar = self.KolmogorovCovar.T * (self.LambdaSeeing/(numpy.pi*2.0))**2.0
        self.OffControl = numpy.real(self.OffControl) * (self.LambdaSeeing/(numpy.pi*2.0))**2.0

    def zernikePropagation(self):
        self.ZernikePropagation = self.Voltage2Zernike[:60,:].T.dot(self.CM[:60,:]).dot(self.HOIM.dot(self.Z2DM))
        self.ZernikePropagation[0][0] = 1.0
        self.ZernikePropagation[1][1] = 1.0
        self.ZernikePropagation[0][1] = 0.0
        self.ZernikePropagation[1][0] = 0.0
        self.PropagatedKolmogorovCovar = self.ZernikePropagation.dot(self.KolmogorovCovar).dot(self.ZernikePropagation.T)

    def seeingEstimation(self, ax=None):
        self.LocalNoiseFactor = numpy.abs(self.RTC/(1.0+self.WFS*self.RTC))**2.0
        self.LocalAtmosphereFactor = numpy.abs((1.0+self.WFS*self.RTC)/(self.WFS*self.RTC))**2.0
        self.A=self.ZPowerCommands.T-numpy.array([self.ZPowerNoise[n,:] * self.LocalNoiseFactor[n] for n 
                     in range(self.LocalNoiseFactor.shape[0])])
        self.ZPowerAtmosphere = numpy.array([self.A[n,:]*self.LocalAtmosphereFactor[n] for n
                     in range(self.A.shape[0])])
        AtmosphereTotalPower = self.ZPowerdFrequencies*numpy.sum(self.ZPowerAtmosphere[:,self.ZIn])
        K = numpy.sum(numpy.array(self.PropagatedKolmogorovCovar.diagonal())[self.ZIn])
        self.SeeingScale = AtmosphereTotalPower/K
        self.Seeing = numpy.real(self.SeeingScale**(3.0/5.0)/self.Arcsec)
        
    def computeSpectralSlope(self, FLim):
        IPropagated = numpy.sum(numpy.array(self.PropagatedKolmogorovCovar.diagonal())[self.ZIn])
        IFull = numpy.sum(numpy.array(self.KolmogorovCovar.diagonal()))
        Ratio = IPropagated/IFull

        Wavelength = 0.5e-6
        TargetIntegral = Ratio * (Wavelength /(2.0*numpy.pi))**2.0

        A = numpy.sum(self.ZPowerCommands[self.ZIn, :], axis=0)
        
        In = (self.ZPowerFrequencies > FLim[0]) & (self.ZPowerFrequencies < FLim[1])

        F = self.ZPowerFrequencies[In]

        Power = A[In]

        LogF = numpy.log10(F)
        LogPower = numpy.log10(Power)

        m = numpy.array([numpy.ones(LogF.shape), LogF])
        im = numpy.linalg.pinv(m)

        self.SpectralSlopes = im.T.dot(LogPower)

        F = self.ZPowerFrequencies
        LogF = numpy.log10(F)
        m = numpy.array([numpy.ones(LogF.shape), LogF])
        Model = 10.0**(m.T.dot(self.SpectralSlopes))
        LowFreq = F <= FLim[0]
        Model[LowFreq] = A[LowFreq]

        dF = numpy.mean(numpy.diff(F))
        dt = 0
        OK = True
        Tau0Max = 0.02
        Step0 = 0.001
        dFModel = Model*dF

        while True:
            dt = dt+ 0.001
            if dt > Tau0Max:
                OK = False
                break
            I = ComputeIntegral(dFModel, F, dt)
            if I > TargetIntegral:
                break
        if not(OK):
            self.Tau0 = Tau0Max
            return
        Low = dt - Step0
        High = dt
        Step = Step0
        while Step >= Step0/20:
            Middle = (Low + High)/2.0
            I = ComputeIntegral(dFModel, F, Middle)
            if I > TargetIntegral:
                High = Middle
            else:
                Low = Middle
            Step = Step / 2.0
        ILow = ComputeIntegral(dFModel, F, Low)
        IHigh = ComputeIntegral(dFModel, F, High)

        interp = scipy.interpolate.interp1d([ILow, IHigh], [Low, High])
        self.Tau0 = interp(numpy.real(TargetIntegral)).tolist()

        return 

    def estimateStrehlRatio(self, ax=None):
        self.OffControl = 1.0 * self.OffControl * self.SeeingScale
        self.ScaledKolmogorovCovar = self.KolmogorovCovar * self.SeeingScale
        self.PropagatedKolmogorovCovar *= self.SeeingScale
        self.UnPropagated = (self.ScaledKolmogorovCovar - self.PropagatedKolmogorovCovar).diagonal()
        self.OffControl += numpy.sum(self.UnPropagated) 

        self.ZPowerWFS = self.ZPowerSlopes.T - self.ZPowerNoise
        self.ZPowerWFS = numpy.array([self.ZPowerWFS[n,:]/(numpy.abs(self.WFS)**2.0)[n] for n in range(self.ZPowerWFS.shape[0])])

        self.TotalWFE2 = numpy.max([0.0, self.WFS2+self.OffControl])
        self.WFE = numpy.sqrt(self.TotalWFE2)

        self.WFSError = numpy.sqrt(numpy.max([0.0, self.WFS2]))
        self.Strehl = numpy.real(numpy.exp(-(2.0*numpy.pi*self.WFE/self.LambdaStrehl)**2.0))
        self.TemporalError = numpy.exp(-(2.0*numpy.pi*self.WFSError/self.LambdaStrehl)**2.0)
        self.rms = numpy.mean(numpy.std(self.Gradients))
        if ax != None:
            print('OffControl: %E' % self.OffControl)
            print('WFE: %E' % self.WFE)
            print('WFSError: %E' % self.WFSError)
            print('Strehl: %f' % self.Strehl)
            raw_input()

    def noiseEvaluation(self):
        FMax = numpy.max(self.ZPowerFrequencies)
        Weights = self.AORejection * self.ZPowerFrequencies/FMax
        Weights[self.ZPowerFrequencies < FMax/2.0] = 0
        Noise = numpy.sum(Weights * self.ZPowerSlopes, axis=1)/numpy.sum(Weights)

        Reference = numpy.array([self.AORejection]).T.dot(numpy.array([numpy.sum(self.S2Z**2.0, axis=1)/FMax]))
        ReferenceNoise = Weights.dot(Reference)/numpy.sum(Weights)

        self.WFSNoise = numpy.sqrt(numpy.mean(Noise[self.ZIn])/numpy.mean(ReferenceNoise[self.ZIn]))
        #self.ZPowerNoise = Reference * self.WFSNoise**2.0
        self.ZPowerNoise = self.WFSNoise**2.0 * numpy.ones([self.ZPowerFrequencies.shape[0],1]).dot(numpy.array([numpy.sum(self.S2Z**2.0, axis=1)]))/numpy.max(self.ZPowerFrequencies)

    def computeResiduals(self):
        Residuals = []
        for frame in self.Gradients:
            Residuals.append(self.S2Z.T.dot(frame))
        self.Residuals = numpy.array(Residuals)

    def computeTTResiduals(self, saveData=False):
        TTR = []
        for frame in self.Gradients:
            tip = numpy.mean(frame[0::2])
            tilt = numpy.mean(frame[1::2])
            TTR.append([tip, tilt])
        self.TTResiduals = numpy.array(TTR)

        if saveData:
            TTResid = numpy.std(self.TTResiduals, axis=0)*0.5
            sqlCommand = "UPDATE CIAO_%d_DataLoggers SET TIP_RESIDUALS = %.4f, TILT_RESIDUALS = %.4f WHERE TIMESTAMP = %.7f;" % (self.CIAO_ID, TTResid[0], TTResid[1], self.startTime)
            self.sqlCursor.execute(sqlCommand)

class CircularBuffer( object ):
    def __init__(self, df='', CDMS_BaseDir = '', CDMS_ConfigDir='', S2M = None, ModalBasis = None, Z2DM = None, S2Z = None, HOIM = None, CM=None, TT2HO=None,
                 DM2Z=None, TTM2Z=None, loopRate=500.0, RTC_Delay=0.5e-3, prefix=''):
        self.df = df
        self.CDMS_BaseDir = CDMS_BaseDir
        self.CDMS_ConfigDir = CDMS_ConfigDir
        self.directory = os.path.dirname(self.df)
        self.header = pyfits.getheader(df)
        self.data = pyfits.getdata(df)
        self.columns = self.data.columns
        self.Intensities = self.data.field('Intensities')
        self.Gradients = self.data.field('Gradients')
        self.HODM = self.data.field('HODM_Positions')
        self.TTM = self.data.field('ITTM_Positions')
        self.FrameCounter = self.data.field('FrameCounter')
        self.time = self.data.field('Seconds')+self.data.field('USeconds')/100000.0
        self.time -= self.time[0]
        try:
            self.loopRate = self.header.get('ESO AOS LOOP RATE')
        except:
            self.loopRate = loopRate
        self.RTC_Delay = RTC_Delay
        self.controlGain = self.header.get('ESO AOS GLOBAL GAIN')
        self.LambdaSeeing = 0.5e-6
        self.LambdaStrehl = 2.2e-6
        self.ReferenceSeeing = 1.0/(180.0*3600.0/numpy.pi)
        self.L0 = numpy.inf
        self.ApertureDiameter = 8.0
        self.r0 = self.LambdaSeeing/self.ReferenceSeeing
        self.CIAO_ID = self.header.get("ESO OCS SYS ID")
        self.Alt = self.header.get("ESO TEL ALT")
        self.Az = self.header.get("ESO TEL AZ")
        if S2M != None:
            self.S2M = pyfits.getdata(S2M)
        else:
            self.S2M = None
        if ModalBasis != None:
            self.ModalBasis = pyfits.getdata(self.CDMS_BaseDir+str(self.CIAO_ID)+self.CDMS_ConfigDir+
                                ModalBasis, ignore_missing_end=True)
        else:
            self.ModalBasis = None
        if HOIM != None:
            self.HOIM = pyfits.getdata(HOIM, ignore_missing_end=True)
        else:
            self.HOIM = None
        if TT2HO != None:
            self.TT2HO = pyfits.getdata(self.CDMS_BaseDir+str(self.CIAO_ID)+self.CDMS_ConfigDir+TT2HO, 
                                        ignore_missing_end=True)
            self.iTT2HO = pyfits.getdata(self.CDMS_BaseDir+str(self.CIAO_ID)+self.CDMS_ConfigDir+
                                        'RecnOptimiser.ITT2HO.fits', ignore_missing_end=True)
        else:
            self.TT2HO = None
            self.iTT2HO = None
        if CM != None:
            self.CM = pyfits.getdata(CM, ignore_missing_end=True)
            self.CM /= self.controlGain
            self.CM[:60,:] += self.TT2HO.dot(self.CM[60:,:])
        else:
            self.CM = None
        if TTM2Z != None:
            self.TTM2Z = pyfits.getdata(self.CDMS_BaseDir+prefix+self.CDMS_ConfigDir+TTM2Z,
                                        ignore_missing_end = True)
        else:
            self.TTM2Z = None
        if Z2DM != None:
            self.Z2DM = pyfits.getdata(self.CDMS_BaseDir+str(self.CIAO_ID)+self.CDMS_ConfigDir+Z2DM,
                                       ignore_missing_end=True)
            self.NZernikes = self.Z2DM.shape[1]
            self.ZIn = numpy.array([(i >= 3) & (i < 99) for i in range(self.NZernikes)])
            if DM2Z == None:
                self.DM2Z = linalg.pinv(self.Z2DM)
            else:
                self.DM2Z = pyfits.getdata(self.CDMS_BaseDir+str(self.CIAO_ID)+self.CDMS_ConfigDir+DM2Z,
                                           ignore_missing_end=True)
        else:
            if DM2Z != None:
                self.DM2Z = pyfits.getdata(self.CDMS_BaseDir+str(self.CIAO_ID)+self.CDMS_ConfigDir+DM2Z, 
                                           ignore_missing_end=True)
                self.Z2DM = linalg.pinv(self.DM2Z)
            try:
                self.Z2DM = pyfits.getdata(os.path.dirname(os.path.realpath(__file__))+
                         '/../data/cimdatZernikeZ2DM.fits', ignore_missing_end=True)
                self.DM2Z = linalg.pinv(self.Z2DM)
            except:
                self.Z2DM = None
                self.DM2Z = None
        if S2Z != None:
            self.S2Z = pyfits.getdata(Z2DM, ignore_missing_end=True)
            self.Z2S = linalg.pinv(self.S2Z)
        else:
            try:
                self.S2Z = self.DM2Z.dot(self.CM[:60])
                self.Z2S = linalg.pinv(self.S2Z)
            except:
                self.S2Z = None
                self.Z2S = None
                raise Exception("inverting CM failed!")
        self.Voltage2Zernike = numpy.concatenate((self.DM2Z.T, self.TTM2Z.T))

    def loadTemplateMaps(self):
        self.expno = self.header.get("HIERARCH ESO TPL EXPNO")
        self.HO_ACT_REF_MAP = pyfits.getdata(self.directory+
                '/HOCtr.ACT_POS_REF_MAP_%04d.fits' %(self.expno+1))
        self.CM = pyfits.getdata(self.directory+
                '/Recn.REC1.CM_%04d.fits' % (self.expno+1))
        self.S2M = pyfits.getdata(self.directory+
                '/RecnOptimiser.S2M_%04d.fits' % (self.expno+1))
        self.M2V = pyfits.getdata(self.directory+
                '/RecnOptimiser.M2V_%04d.fits' % (self.expno+1))
        self.V2M = linalg.pinv(self.M2V)
        self.getDisturbanceRegions()

    def getDisturbanceRegions(self, start = None, stop = None):
        if start == None:
            disturbance = []
            calm = []
            for i in range(len(self.HODM)):
                if (self.HODM[i] == self.HO_ACT_REF_MAP).all():
                    calm.append(i)
                else:
                    disturbance.append(i)
            self.disturbance = numpy.array(disturbance)
            self.calm = numpy.array(calm)
        else:
            disturbance = numpy.array(((self.FrameCounter > start) & (self.FrameCounter \
                    <= stop)))
            self.disturbance = numpy.arange(len(disturbance))[disturbance]
            self.calm = numpy.arange(len(disturbance))[disturbance==False]


    def getAberrations(self):
        calmSlopes = numpy.average(self.Gradients[self.calm], axis=0)
        meanActPos = self.CM.dot(calmSlopes)
        #self.aberrations = self.calculateMirrorZernikes(self.HO_ACT_REF_MAP[0]
        #        + meanActPos[:-2])
        self.aberrations = self.calculateMirrorZernikes(meanActPos[:-2])
        return self.aberrations
        
        

    def calculateMirrorZernikes(self, MirrorPosition):
        self.MirrorZernikes = self.DM2Z.T.dot(MirrorPosition)
        return self.MirrorZernikes

    def calculatePSD(self):
        self.modes = []
        for frame in self.Gradients:
            self.modes.append(self.S2M.dot(frame))

        self.modes = numpy.array(self.modes).T
        self.FFT = fftpack.fftshift(fftpack.fft(self.modes))
        self.FFT = self.FFT * self.FFT.conjugate()
        self.freq = fftpack.fftshift(fftpack.fftfreq(len(self.modes[0]), 1.0/526.7))
    
    def processFOV(self, ax=None):
        FOV = FOVFrame()
        startingPoint = self.TTM[0]
        for tt, frame in zip(self.TTM, self.Intensities):
            if ((tt[0] != startingPoint[0]) & (tt[1] != startingPoint[1])):
                FOV.addFrame(tt[0], tt[1], frame)

        FOV.publish(ax)

    def extractModulation(self, ax = None):
        firstActuator = self.HODM[:,0]
        modulation = numpy.arange(len(firstActuator))[firstActuator != firstActuator[0]]
        self.modulation_start = modulation[0]
        self.modulation_stop = modulation[-1]
        self.modes = []
        if ax != None:
            ax.clear()
        for frame in self.HODM[self.modulation_start:self.modulation_stop]:
            self.modes.append(self.ModalBasis.dot(frame))
            if ax != None:
                ax.plot(self.modes[-1])
        self.modes = numpy.array(self.modes)
        if ax != None:
            ax.figure.show()
            raw_input('Press Enter to Continue')

    def extractIM(self, ax = None):
        counter = 0
        IM = []
        for m in self.HODM[self.modulation_start:self.modulation_stop].T:
            response = []
            counter += 1
            for slope in self.Gradients[self.modulation_start+1:self.modulation_stop+1,:].T:
                response.append(signal.correlate(slope, m, mode='valid')[0])
            IM.append(response)
        
        self.IM = numpy.array(IM)
        if ax != None:
            ax.clear()
            ax.imshow(self.IM.T)
            ax.figure.show()
            raw_input()

    def processSlopes(self):
        self.ProcessedGradients = ProcessedData(self.Gradients, self.S2M, self.S2Z)
        
    def processVoltages(self):
        self.ProcessedVoltages = ProcessedData(self.HODM, self.V2M[:,:-2]/1e-6,
                self.DM2Z.T/1e-6)
        
    def synchronizeData(self):
        self.TTM = self.TTM[:-3,:]
        #self.TTM = self.TTM[2:-1,:]
        self.HODM = self.HODM[2:-1,:]
        self.Gradients = self.Gradients[3:,:]
        
        self.TTM -= numpy.mean(self.TTM, axis=0)
        self.HODM -= numpy.mean(self.HODM, axis=0)
        self.Gradients -= numpy.mean(self.Gradients, axis=0)
        
    def zernikeSpace(self):
        self.ZSlopes = self.Gradients.dot(self.S2Z.T)
        self.ZCommands = numpy.concatenate((self.HODM.T, self.TTM.T)).T.dot(self.Voltage2Zernike)
        
    def computePSD(self, source):
        NWindows = 1
        Resolution = 0.5
        if source == 'ZSlopes':
            data = self.ZSlopes
        elif source == 'ZCommands':
            data = self.ZCommands
        NWindows = computeNWindows(data, self.loopRate, Resolution)
        NSeries = data.shape[1]
        data = FourierDetrend(data, 'TwoPoints') 

        NFrames = data.shape[0]
        Power = numpy.abs(numpy.fft.fft(data.T/NFrames))**2.0
        Freq = numpy.fft.fftfreq(NFrames, d=1.0/self.loopRate)

        Power = numpy.fft.fftshift(Power, axes=1)
        Freq = numpy.fft.fftshift(Freq)
        Power = 2.0*Power[:,Freq >= 0]
        Freq = Freq[Freq >= 0]
        Power[:,0] /= 2.0
        Power = Power / numpy.mean(numpy.diff(Freq))

        NSamplesPerWindow = numpy.floor(Power.shape[1]/NWindows)
        TotalSampleNumber = NWindows * NSamplesPerWindow
        Power = Power.T

        if NWindows > 1:
            Power = Power[:TotalSampleNumber, :]
            Power = Power.reshape([NWindows, NSamplesPerWindow, NSeries], order='F')
            Power = numpy.mean(Power, axis=0)

            Freq = Freq[:TotalSampleNumber]
            Freq = Freq.reshape([NWindows, NSamplesPerWindow], order = 'F')
            Freq = numpy.mean(Freq, axis = 0)

        if source == 'ZSlopes':
            self.ZPowerSlopes = Power
            self.ZPowerFrequencies = Freq
            self.ZPowerdFrequencies = numpy.mean(numpy.diff(self.ZPowerFrequencies))
        elif source == 'ZCommands':
            self.ZPowerCommands = Power
        return

    def AORejectionFunction(self):     # AORejection, WFS, RTC, NoisePropagation
        self.WFS2 = numpy.sum(self.ZPowerSlopes) * self.ZPowerdFrequencies
        Omega = 2*numpy.pi*self.ZPowerFrequencies
        iOmega = (0.0 + 1.0j)*Omega
        
        RTC_Period = 1.0/self.loopRate
        CycleDelay = numpy.exp(-RTC_Period*iOmega)
        CPU_Delay = numpy.exp(-self.RTC_Delay*iOmega)

        GainNumerator = numpy.array([self.controlGain])
        GainDenominator = numpy.array([1, -1])
        Numerator = GetTF(GainNumerator, CycleDelay)
        Denominator = GetTF(GainDenominator, CycleDelay)
        Controller = Numerator/Denominator


        DM = 1

        Ratio = self.RTC_Delay/RTC_Period
        N = numpy.floor(Ratio)
        alpha = Ratio - N

        RSinPhi = -alpha*numpy.sin(Omega*RTC_Period)
        RCosPhi = 1.0-alpha*(1.0-numpy.cos(Omega*RTC_Period))
        Phi = numpy.arctan2(RSinPhi, RCosPhi)
        R = numpy.sqrt(RSinPhi**2.0 + RCosPhi**2.0)

        self.WFS = R*numpy.exp((0.0+1.0j)*Phi)*CycleDelay**(N+1)

        self.RTC = Controller * DM
        OpenLoop_TF = self.WFS*self.RTC
        self.AORejection = 1.0/(1.0+OpenLoop_TF)
        self.NoisePropagation = self.RTC*self.AORejection

        self.AORejection[numpy.isnan(self.AORejection)]= 0.0
        self.WFS[numpy.isnan(self.WFS)] = 0.0
        self.RTC[numpy.isnan(self.RTC)] = 0.0
        self.NoisePropagation[numpy.isnan(self.NoisePropagation)] = 0.0
        self.AORejection = numpy.abs(self.AORejection)**2.0

    #def measureVibrations(self)

    def combinePSDs(self):
        RTC2 = numpy.abs(self.RTC)**2.0

        RTC2S = numpy.array([RTC2 * self.ZPowerSlopes[:,n] for n in range(self.ZPowerSlopes.shape[1])]).T
        self.ZPowerCommands = (RTC2S.T + [RTC2*self.ZPowerCommands[:,n] for n in range(self.ZPowerCommands.shape[1])])/(1.0+RTC2).T
        self.ZPowerSlopes = self.ZPowerCommands/RTC2
        

    def computeKolmogorovCovar(self):
        D_L0 = self.ApertureDiameter/self.L0

        infOuterScaleCovar, radial_orders = newKTaiaj(self.NZernikes+1, self.NZernikes+1, self.ApertureDiameter, self.r0)
        infOuterScaleCovar[:,0] = 0.0
        infOuterScaleCovar[0,:] = 0.0

        max_radial_order = numpy.int(numpy.max(radial_orders))

        reduc_var = OuterScale2VarianceReduction(D_L0, max_radial_order)
        reduction_factors = numpy.ones(radial_orders.shape)

        interp = scipy.interpolate.interp1d(numpy.array(range(max_radial_order))+1, reduc_var)
        blah = radial_orders != 0.0
        reduction_factors[blah] = interp(radial_orders[blah])

        self.KolmogorovCovar = infOuterScaleCovar * reduction_factors
        FullVariance, n = KolmogorovVariance(self.ApertureDiameter, self.L0, self.r0, 1000)

        self.FullVariance = numpy.sum(FullVariance)
        self.FractionOfTotalVariance = numpy.real(numpy.sum(numpy.diag(self.KolmogorovCovar))/self.FullVariance)

        self.OffControl = numpy.sum(numpy.diag(self.KolmogorovCovar))*(1.0/self.FractionOfTotalVariance -1.0)
        self.KolmogorovCovar = self.KolmogorovCovar[1:, 1:]

        self.KolmogorovCovar = self.KolmogorovCovar.T * (self.LambdaSeeing/(numpy.pi*2.0))**2.0
        self.OffControl = numpy.real(self.OffControl) * (self.LambdaSeeing/(numpy.pi*2.0))**2.0

    def zernikePropagation(self):
        self.ZernikePropagation = self.Voltage2Zernike[:60,:].T.dot(self.CM[:60,:]).dot(self.HOIM.dot(self.Z2DM))
        self.ZernikePropagation[0][0] = 1.0
        self.ZernikePropagation[1][1] = 1.0
        self.ZernikePropagation[0][1] = 0.0
        self.ZernikePropagation[1][0] = 0.0
        self.PropagatedKolmogorovCovar = self.ZernikePropagation.dot(self.KolmogorovCovar).dot(self.ZernikePropagation.T)

    def seeingEstimation(self, ax=None):
        self.Arcsec = 180.0*3600.0/numpy.pi
        self.LocalNoiseFactor = numpy.abs(self.RTC/(1.0+self.WFS*self.RTC))**2.0
        self.LocalAtmosphereFactor = numpy.abs((1.0+self.WFS*self.RTC)/(self.WFS*self.RTC))**2.0
        self.A=self.ZPowerCommands.T-numpy.array([self.ZPowerNoise[n,:] * self.LocalNoiseFactor[n] for n 
                     in range(self.LocalNoiseFactor.shape[0])])
        self.ZPowerAtmosphere = numpy.array([self.A[n,:]*self.LocalAtmosphereFactor[n] for n
                     in range(self.A.shape[0])])
        AtmosphereTotalPower = self.ZPowerdFrequencies*numpy.sum(self.ZPowerAtmosphere[:,self.ZIn])
        K = numpy.sum(numpy.array(self.PropagatedKolmogorovCovar.diagonal())[self.ZIn])
        self.SeeingScale = AtmosphereTotalPower/K
        self.Seeing = numpy.real(self.SeeingScale**(3.0/5.0)/self.Arcsec)
        
    def computeSpectralSlope(self, FLim):
        IPropagated = numpy.sum(numpy.array(self.PropagatedKolmogorovCovar.diagonal())[self.ZIn])
        IFull = numpy.sum(numpy.array(self.KolmogorovCovar.diagonal()))
        Ratio = IPropagated/IFull

        Wavelength = 0.5e-6
        TargetIntegral = Ratio * (Wavelength /(2.0*numpy.pi))**2.0

        A = numpy.sum(self.ZPowerCommands[self.ZIn, :], axis=0)
        
        In = (self.ZPowerFrequencies > FLim[0]) & (self.ZPowerFrequencies < FLim[1])

        F = self.ZPowerFrequencies[In]

        Power = A[In]

        LogF = numpy.log10(F)
        LogPower = numpy.log10(Power)

        m = numpy.array([numpy.ones(LogF.shape), LogF])
        im = numpy.linalg.pinv(m)

        self.SpectralSlopes = im.T.dot(LogPower)

        F = self.ZPowerFrequencies
        LogF = numpy.log10(F)
        m = numpy.array([numpy.ones(LogF.shape), LogF])
        Model = 10.0**(m.T.dot(self.SpectralSlopes))
        LowFreq = F <= FLim[0]
        Model[LowFreq] = A[LowFreq]

        dF = numpy.mean(numpy.diff(F))
        dt = 0
        OK = True
        Tau0Max = 0.02
        Step0 = 0.001
        dFModel = Model*dF

        while True:
            dt = dt+ 0.001
            if dt > Tau0Max:
                OK = False
                break
            I = ComputeIntegral(dFModel, F, dt)
            if I > TargetIntegral:
                break
        if not(OK):
            self.Tau0 = Tau0Max
            return
        Low = dt - Step0
        High = dt
        Step = Step0
        while Step >= Step0/20:
            Middle = (Low + High)/2.0
            I = ComputeIntegral(dFModel, F, Middle)
            if I > TargetIntegral:
                High = Middle
            else:
                Low = Middle
            Step = Step / 2.0
        ILow = ComputeIntegral(dFModel, F, Low)
        IHigh = ComputeIntegral(dFModel, F, High)

        interp = scipy.interpolate.interp1d([ILow, IHigh], [Low, High])
        self.Tau0 = interp(numpy.real(TargetIntegral)).tolist()

        return 

    def estimateStrehlRatio(self, ax=None):
        self.OffControl = 1.0 * self.OffControl * self.SeeingScale
        self.ScaledKolmogorovCovar = self.KolmogorovCovar * self.SeeingScale
        self.PropagatedKolmogorovCovar *= self.SeeingScale
        self.UnPropagated = (self.ScaledKolmogorovCovar - self.PropagatedKolmogorovCovar).diagonal()
        self.OffControl += numpy.sum(self.UnPropagated) 

        self.ZPowerWFS = self.ZPowerSlopes.T - self.ZPowerNoise
        self.ZPowerWFS = numpy.array([self.ZPowerWFS[n,:]/(numpy.abs(self.WFS)**2.0)[n] for n in range(self.ZPowerWFS.shape[0])])

        self.TotalWFE2 = numpy.max([0.0, self.WFS2+self.OffControl])
        self.WFE = numpy.sqrt(self.TotalWFE2)

        self.WFSError = numpy.sqrt(numpy.max([0.0, self.WFS2]))
        self.Strehl = numpy.real(numpy.exp(-(2.0*numpy.pi*self.WFE/self.LambdaStrehl)**2.0))
        self.TemporalError = numpy.exp(-(2.0*numpy.pi*self.WFSError/self.LambdaStrehl)**2.0)
        self.rms = numpy.mean(numpy.std(self.Gradients))
        if ax != None:
            print('OffControl: %E' % self.OffControl)
            print('WFE: %E' % self.WFE)
            print('WFSError: %E' % self.WFSError)
            print('Strehl: %f' % self.Strehl)
            raw_input()

    def noiseEvaluation(self):
        FMax = numpy.max(self.ZPowerFrequencies)
        Weights = self.AORejection * self.ZPowerFrequencies/FMax
        Weights[self.ZPowerFrequencies < FMax/2.0] = 0
        Noise = numpy.sum(Weights * self.ZPowerSlopes, axis=1)/numpy.sum(Weights)

        Reference = numpy.array([self.AORejection]).T.dot(numpy.array([numpy.sum(self.S2Z**2.0, axis=1)/FMax]))
        ReferenceNoise = Weights.dot(Reference)/numpy.sum(Weights)

        self.WFSNoise = numpy.sqrt(numpy.mean(Noise[self.ZIn])/numpy.mean(ReferenceNoise[self.ZIn]))
        #self.ZPowerNoise = Reference * self.WFSNoise**2.0
        self.ZPowerNoise = self.WFSNoise**2.0 * numpy.ones([self.ZPowerFrequencies.shape[0],1]).dot(numpy.array([numpy.sum(self.S2Z**2.0, axis=1)]))/numpy.max(self.ZPowerFrequencies)

    def computeTTResiduals(self):
        TTR = []
        for frame in self.Gradients:
            tip = numpy.mean(frame[0::2])
            tilt = numpy.mean(frame[1::2])
            TTR.append([tip, tilt])
        self.TTResiduals = numpy.array(TTR)

def ComputeIntegral(dFModel, F, dt):
    Integrant = dFModel * numpy.abs(1.0 - numpy.cos(2.0*numpy.pi*F*dt))
    return numpy.sum(Integrant)

def KolmogorovVariance(D, L0, r0, nOrders):
    cst = gamma(11./6.)**2.0*gamma(14./3.)*(24.0*gamma(6./5.)/5.)**(5./6.)/(2.**(8./3.)*numpy.pi)/gamma(17./6.)**2.

    n, m, Sign = FindNM(nOrders)

    y = cst * (-1)**(n-m)*(n+1.0)*gamma(n-5./6.)/gamma(n+23./6.)

    reduc_var = OuterScale2VarianceReduction(D/L0,numpy.max(n))
    reduction_factors = numpy.ones(n.shape)
    blah = n != 0.0
    for i in range(nOrders):
        if blah[i]:
            reduction_factors[i] = reduc_var[n[i]-1]
    #reduction_factors[blah] = reduc_var[n[blah]]

    y = y * reduction_factors
    y = y * (D/r0) **(5./3.)
    y[0] = 0
    i = numpy.array(range(nOrders))
    y = y[i]

    return y, n

    

def OuterScale2VarianceReduction(D, maxrad):
    n_rad = numpy.array(range(int(maxrad)))+1.0

    x = D/2.
    y = x**(2*n_rad-5.0/3.0)
    w = 1

    p = 0
    facto_p = 1
    sign_p = 1
    reduc_var = 0

    mult = gamma(n_rad-5./6.) / gamma(n_rad + 23./6.)
    variance = 0.756 * (n_rad+1)*mult

    increase = 1
    threshold = 1e-5

    while numpy.max(numpy.abs(increase)) > threshold:
        mult = gamma(p+n_rad+1.5)*gamma(-p-n_rad+5./6.)*gamma(p+n_rad+1.0)
        mult /= gamma(p+2*n_rad+3.0)*gamma(p+n_rad+2.0)
        increase = sign_p / facto_p * y * mult

        mult = gamma(-p+n_rad-5.0/6.0) * gamma(p+7.0/3.0) * gamma(p+11.0/6.0)
        mult /= gamma(p+n_rad+23.0/6.0)*gamma(p+17./6.)
        increase = increase + sign_p /facto_p*w*mult

        reduc_var = reduc_var + increase
        increase = increase / reduc_var
        y = y*x**2.0
        w = w*x**2.0
        p = p+1
        facto_p = facto_p*p
        sign_p = -sign_p

    reduc_var = reduc_var * 1.1641 * (n_rad+1.0)/variance
    reduc_var[n_rad<1] = 1.0
    reduc_var[reduc_var > 1] = 1.0

    return reduc_var

def newKTaiaj(i, j, D, r0):
    cst = gamma(11.0/6.0)**2.0*gamma(14./3.)*(24*gamma(6.0/5.0)/5.0)**(5.0/6.0)/(2.0**(8.0/3.0)*numpy.pi)

    ni, mi, Sign = FindNM(i)
    nj, mj, Sign = FindNM(j)

    NI = []
    MI = []
    NJ = []
    MJ = []
    I = []
    J = []
    for x in range(j):
        NI.append(ni)
        MI.append(mi)
        NJ.append(nj)
        MJ.append(mj)
        I.append(numpy.array(range(i)) + 1.0)
        J.append(numpy.array(range(j)) + 1.0)

    ni = numpy.array(NI)
    mi = numpy.array(MI)
    nj = numpy.array(NJ).T
    mj = numpy.array(MJ).T
    i = numpy.array(I)
    j = numpy.array(J).T

    n = numpy.sqrt(ni*nj)
    y = cst *(-1.0+0j)**((ni+nj-2*mi)/2.0)*numpy.sqrt((ni+1)*(nj+1))*gamma((ni+nj-5./3.)/2.) / (gamma((ni-nj+17./3.)/2.0) * gamma((nj-ni+17.0/3.0)/2.0) *
              gamma((ni+nj+23./3.0)/2.0))
    y[mi!=mj] = 0.0
    blah = numpy.remainder(numpy.abs(i-j),2)==1
    junk = mi !=0
    y *= 1.0*((blah & junk) == False)
    y = y * (D/r0)**(5.0/3.0)

    return y, n

def FindNM(i):
    mode = numpy.array(range(i))+1.0
    nf = numpy.ceil(numpy.sqrt(2.0*mode+0.25)-1.5)

    mf = mode - nf*(nf+1)/2.0
    odd = numpy.mod(nf,2) == 1

    mf[odd] = 2*numpy.floor((mf[odd]+1)/2) -1
    mf[odd==False] = 2*numpy.floor(mf[odd==False]/2)
    Sign = numpy.mod(mode - (nf-1)*nf/2,2)

    return nf, mf, Sign



def GetTF(Coeffs, Delay):
    H = 0
    for k in range(len(Coeffs)):
        H = H + Coeffs[k]*Delay**(k)
    return H

def computeNWindows(data, loopRate, R):
    N = data.shape[0]
    dF = loopRate/(1.0*N)
    NWindows = round(R/dF)
    NWindows = max(1, NWindows)
    return NWindows


def FourierDetrend(data, detrend):
    if detrend.lower() == 'twopoints':
        n = data.shape[0]
        m = numpy.ones((n,2))
        m[:,1] = [j+1-(1+n)/2.0 for j in range(n)]
        im = linalg.pinv(m[(0,-1),:])
        c = im.dot(data[(0,-1),:])

        data = data-m.dot(c)

        return  data - numpy.mean(data, axis=0)


class ProcessedData ( object ):
    def __init__(self, data, modeProjection, zernikeProjection):
        self.average = numpy.average(data, axis=0)
        modes = []
        zerns = []
        for frame in data:
            modes.append(modeProjection.dot(frame))
            zerns.append(zernikeProjection.dot(frame))
        self.modes = numpy.array(modes)
        self.modes_average = numpy.average(self.modes, axis=0)
        self.modes_std = numpy.std(self.modes, axis=0)
        self.zerns = numpy.array(zerns)
        self.zerns_average = numpy.average(self.zerns, axis=0)
        self.zerns_std = numpy.std(self.zerns, axis=0)

class Controller( object ):
    def __init__(self, CM = None):
        self.CM = CM
        self.HOCM = self.CM[:60,:]
        self.TTCM = self.CM[60:,:]
        self.iTT2HO = pyfits.getdata('../../data/RecnOptimiser.ITT2HO.fits',
                ignore_missing_end=True)
        self.TT2HO = pyfits.getdata('../../data/RecnOptimiser.TT2HO.fits',
                ignore_missing_end=True)
        self.SMA = pyfits.getdata('../../data/HOCtr.SMA_MATRIX.fits',
                ignore_missing_end=True)
        self.ModalBasis = pyfits.getdata('../../data/RecnOptimiser.ModalBasis.fits',
                ignore_missing_end=True)
        self.P2DM = pyfits.getdata('../../data/HOCtr.P2DM_PROJECTION.fits',
                ignore_missing_end=True)
        self.iP2DM = pyfits.getdata('../../data/HOCtr.IP2DM_PROJECTION.fits',
                ignore_missing_end=True)
        self.GarbageProj = pyfits.getdata('../../data/HOCtr.GARBAGE_PROJECTION.fits',
                ignore_missing_end=True)
        self.KT = 0.001
        self.KI = 0.5

    def computeDeltas(self, slopes):
        return self.HOCM.dot(slopes)

    def computeCommands(self, slopes):
        #HODM = self.HOCM.dot(slopes) - self.TT2HO.dot(
        #        self.iTT2HO.dot(self.HOCM.dot(slopes)))
        HODM = self.HOCM.dot(slopes)
        TTM = self.TTCM.dot(slopes)
        HODM = HODM + self.TT2HO.dot(TTM)

        return HODM

    def doSMA(self, command):
        saturated = command.copy()
        #saturated[saturated < self.
        return command
        

class LoopAnalyzer( object):
    def __init__(self, df =None):
        self.CB = CircularBuffer(df = df)
        self.CB.loadTemplateMaps()
        self.HOCtr = Controller(self.CB.CM)
        self.predictions = []


    def predict(self):
        predictions = []
        for frame, mirror in zip(self.CB.Gradients, self.CB.HODM):
            Pipeline = (1.0-self.HOCtr.KT)*mirror - \
                      self.HOCtr.KI * self.HOCtr.computeCommands(frame)
            Pipeline = Pipeline + self.CB.HO_ACT_REF_MAP
            Pipeline = self.HOCtr.doSMA(Pipeline)
            predictions.append(Pipeline)
        
        self.predictions = numpy.array(predictions)

class AcqCamImage( object ):
    def __init__(self, datadir='', df=''):
        self.datadir = datadir
        self.df = df
        self.imageCube = pyfits.getdata(self.datadir+self.df)
        self.sky = pyfits.getdata(self.datadir+'sky.fits')
        self.dead = pyfits.getdata(self.datadir+'dead.fits')

        self.x = numpy.array([1565, 1130, 660, 150])
        self.y = numpy.array([130, 120, 125, 125])
        self.width = 250
        self.clean()
        #self.mean = numpy.mean(self.imageCube)
        #self.std = numpy.std(self.imageCube)

    def clean(self):
        self.cleanSky = {}
        for ut in range(4):
            cleaned = numpy.zeros(self.sky.shape)
            for x in numpy.arange(self.x[ut], self.x[ut]+self.width):
                for y in numpy.arange(self.y[ut], self.y[ut]+self.width):
                    if self.dead[y][x] == 1:
                        r = 1
                        while True:
                            region = self.dead[y-(r+1):y+r, x-(r+1):x+r]
                            if numpy.sum(region) - region.shape[0]*region.shape[1] < -3.0:
                                break
                            else:
                                r += 1
                        cleaned[y][x] = numpy.median(self.sky[y-(r+1):y+r, x-(r+1):x+r][self.dead[y-(r+1):y+r,x-(r+1):x+r]==0])
                    else:
                        cleaned[y][x] = self.sky[y][x]
            self.cleanSky[ut] = cleaned[self.y[ut]:self.y[ut]+self.width, self.x[ut]:self.x[ut]+self.width].copy()
        #"""
        self.clean = {}
        self.mean = {}
        self.std = {}
        for ut in range(4):
            print("%d" % ut)
            cleaned = []
            mean = []
            std = []
            for frame in self.imageCube:
                for x in numpy.arange(self.x[ut], self.x[ut]+self.width):
                    for y in numpy.arange(self.y[ut], self.y[ut]+self.width):
                        if self.dead[y][x] == 1:
                            r = 1
                            while True:
                                region = self.dead[y-(r+1):y+r, x-(r+1):x+r]
                                if numpy.sum(region) - region.shape[0]*region.shape[1] < -3.0:
                                    break
                                else:
                                    r += 1
                            frame[y][x] = numpy.median(frame[y-(r+1):y+r,x-(r+1):x+r][self.dead[y-(r+1):y+r, x-(r+1):x+r]==0])
                cleaned.append(frame[self.y[ut]:self.y[ut]+self.width, self.x[ut]:self.x[ut]+self.width]- self.cleanSky[ut])
                mean.append(numpy.mean(cleaned[-1]))
                std.append(numpy.std(cleaned[-1]))

            self.clean[ut] = numpy.array(cleaned)
            self.mean[ut] = numpy.array(mean)
            self.std[ut] = numpy.array(std)
        #"""

    def findPeaks(self, size=4):
        neighborhood = generate_binary_structure(2, 2)
        self.x = {}
        self.y = {}
        for ut in range(4):
            for frame in self.clean[ut]:
                local_max = maximum_filter(frame, size=size)==frame
                coordx, coordy =numpy.meshgrid(numpy.arange(frame.shape[0]),numpy.arange(frame.shape[1]))
            self.x[ut] = coordx[local_max]
            self.y[ut] = coordy[local_max]

    def stack(self):
        self.postageStamps = {}
        self.speckleStamps = {}
        for ut in range(4):
            postageStamp = []
            speckleStamp = []
            SS = {}
            for frame, mean, std in zip(self.clean[ut], self.mean[ut], self.std[ut]):
                PS = []
                goodX = []
                goodY = []
                for x, y in zip(self.x[ut], self.y[ut]):
                    if (frame[y,x] > mean+2.0*std) & (10 < x <240) & (10 < y < 240):
                        PS.append(frame[y-10:y+10, x-10:x+10])
                        PS[-1] /= numpy.max(PS[-1])
                        goodX.append(x)
                        goodY.append(y)
                        if str([x,y]) in SS.keys():
                            SS[str([x,y])].append(frame[y-10:y+10, x-10:x+10])
                        else:
                            SS[str([x,y])] = [frame[y-10:y+10, x-10:x+10]]

                postageStamp.append(numpy.median(numpy.array(PS), axis=0))
            for key in SS.keys():
                SS[key] = numpy.array(SS[key])
            self.postageStamps[ut] = numpy.array(postageStamp)
            self.speckleStamps[ut] = SS


class FLIRCamImage( object ):
    def __init__(self, df):
        self.df = df
        self.image = Image.open(self.df)
        self.imdata = numpy.array(self.image)[:,:,:3].sum(axis=2)
        self.mean = numpy.mean(self.imdata[self.imdata > 0])
        self.max = numpy.max(self.imdata)
        self.imdata = self.imdata/float(self.max)
        #self.imdata /= numpy.mean(self.imdata[self.imdata > 0])
        self.floor = numpy.min(self.imdata[self.imdata > 0])
        self.nx = len(self.imdata[0])
        self.ny = len(self.imdata)
        self.npix = self.nx*self.ny

    def Spot(self, p):
        sig = p[0]
        xc = p[1]
        yc = p[2]

        if (p[0] <= 0.0):
            return numpy.zeros(self.ncutout)
        if (p[1] < 0) or (p[1] > self.cutout_x):
            return numpy.zeros(self.ncutout)
        if (p[2] < 0) or (p[2] > self.cutout_y):
            return numpy.zeros(self.ncutout)
        x = numpy.arange(self.cutout_x)
        y = numpy.arange(self.cutout_y)
        image = numpy.zeros((self.cutout_x, self.cutout_x))

        for i in range(self.cutout_x):
            for j in range(self.cutout_y):
                image[i][j] = 1.0*( numpy.exp(-(x[i]-xc)**2.0/sig)*
                        numpy.exp(-(y[j]-yc)**2.0/sig))

        maxval = numpy.max(image)
        image = image/maxval
        image[image < self.floor] = 0.0
        return image.ravel()

    def findFocusCenter(self, x, y, ax=None):
        """
        cutout should be 8x8 grid of pixels
        """

        cutout = self.extractCutout(x, y, 10)

        if ax!=None:
            ax.matshow(cutout)
            ax.figure.show()
            raw_input()
        self.cutout_x = len(cutout[0])
        self.cutout_y = len(cutout)
        self.ncutout = self.cutout_x*self.cutout_y

        sig = 2.4
        xcenter = self.cutout_x/2.0
        ycenter = self.cutout_y/2.0

        errfunc = lambda p, y : numpy.abs(self.Spot(p) - y)
        coeffs = [sig, xcenter, ycenter]
        pfit, success = optimize.leastsq(errfunc, coeffs, args=(cutout.ravel()))

        xfit = pfit[2] - xcenter + x
        yfit = pfit[1] - ycenter + y
        sigma = pfit[0]**0.5
        return sigma, xfit, yfit, cutout

    def findPupilCenter(self, x, y, zoomFactor, pupilImage=None, ax=None):
        """
        cutout should be 8x8 grid of pixels
        """

        cutout = self.extractCutout(x, y, 80, chopTop=True)

        cutoutZoom = self.zoomIn(cutout, zoomFactor)

        corr = signal.correlate2d(cutoutZoom, pupilImage, boundary='fill', mode='valid')
        xc, yc = numpy.unravel_index(numpy.argmax(corr), corr.shape)
        xzp = corr.shape[0]/2.0
        yzp = corr.shape[1]/2.0
        xcent = (xc-xzp)/zoomFactor + x
        ycent = (yc-yzp)/zoomFactor + y

        if ax!=None:
            ax.clear()
            ax.matshow(corr)
            ax.figure.show()
            raw_input()
        return xcent, ycent

    def extractCutout(self, x, y, size, chopTop=False):
        cutout = self.imdata[y-size:y+size, x-size:x+size].copy()

        if chopTop:
            cutoff = numpy.mean(cutout)
            cutout[cutout > cutoff] = 1.0
        return cutout

    def zoomIn(self, original, factor):
        x = numpy.arange(original.shape[0])
        y = numpy.arange(original.shape[1])
        interpolator = scipy.interpolate.RectBivariateSpline(x, y, original)
        new_x = numpy.linspace(0, x[-1], num=len(x)*factor)
        new_y = numpy.linspace(0, y[-1], num=len(y)*factor)
        xx, yy = numpy.meshgrid(new_x, new_y)
        zoomed = interpolator.ev(xx, yy).reshape(len(x)*factor, len(y)*factor)
        return zoomed

class NGCImage( object):
    def __init__(self, filename):
        self.filename = filename
        self.header = pyfits.getheader(filename)
        self.data = pyfits.getdata(filename)

    def subtract_background(self, background):
        self.subtracted = self.data - background.data

    def generate_model(self, n):
        feature_width = int(n*self.spacing+4*self.fwhm)
        feature_x = numpy.arange(feature_width)
        feature_y = numpy.zeros(feature_width)

        for i in range(n):
            c = (i+1)*self.spacing
            feature_y += self.height*numpy.exp(-(feature_x-c)**2.0/(2.0*self.fwhm/2.4)**2.0)

        return feature_x, feature_y

    def find_centers(self, n, collapse):
        distance, height = self.generate_model(n)
        y_corr = scipy.correlate(collapse, height)
        x_corr = scipy.linspace(0, len(y_corr)-1, num =len(y_corr))

        peak = x_corr[numpy.argsort(y_corr)[-1]]

        centers = []
        for i in range(n):
            centers.append((i+1)*self.spacing+peak)

        return numpy.array(centers)

    def findSubapertureCenters(self, nx=20, ny=21):
        self.nx = nx
        self.ny = ny
        self.xcollapse = self.subtracted.sum(axis=0)
        self.ycollapse = self.subtracted.sum(axis=1)

        self.fwhm = 2.5
        self.spacing = 8.1
        self.height = 0.75*numpy.max(self.xcollapse)

        self.xcenters = self.find_centers(nx, self.xcollapse)
        self.ycenters = self.find_centers(ny, self.ycollapse)

    def TopHat(self, p):
        sig_x = p[0]
        sig_y = p[1]
        theta = p[2]
        amplitude = numpy.deg2rad(p[3])
        xc = p[4]
        yc = p[5]

        if (sig_x < 0):# or (sig_x > 18.0):
            return numpy.zeros(64)
        if (sig_y < 0):# or (sig_y > 18.0):
            return numpy.zeros(64)
        fine_grid = numpy.linspace(0, 8)
        fine_image = numpy.zeros((len(fine_grid), len(fine_grid)))

        a = numpy.cos(theta)**2/(2.0*sig_x**2.0) + numpy.sin(theta)**2./(2.0*sig_y**2.0)
        b = -numpy.sin(2.0*theta)**2/(4.0*sig_x**2.0) + numpy.sin(2.0*theta)**2./(4.0*sig_y**2.0)
        c = numpy.sin(theta)**2/(2.0*sig_x**2.0) + numpy.cos(theta)**2./(2.0*sig_y**2.0)

        indices = range(len(fine_grid))
        
        for x in indices:
            for y in indices:
                fine_image[x][y] = amplitude*numpy.exp(-(a*(fine_grid[x]-xc)**2.0 +
                    2.0*b*(fine_grid[x]-xc)*(fine_grid[y]-yc) +
                    c*(fine_grid[y]-yc)**2.0))

        coarse_image = numpy.zeros((8,8))
        for x in range(8):
            for y in range(8):
                goodx = numpy.arange(len(fine_grid))[(x <= fine_grid) & (fine_grid < x+1)]
                goody = numpy.arange(len(fine_grid))[(y <= fine_grid) & (fine_grid < y+1)]
                coarse_image[x][y] = numpy.sum(fine_image[goodx][:,goody])

        return coarse_image.ravel()

    def Elipse(self, p):
        sig_x = p[0]
        sig_y = sig_x*p[1]
        theta = numpy.deg2rad(p[2])
        xc = p[3]
        yc = p[4]

        if (p[1] >= 1.0):
            return numpy.zeros(64)
        if (p[2] < -180.0) or (p[2] > 180.0):
            return numpy.zeros(64)
        if (sig_x < 0):# or (sig_x > 18.0):
            return numpy.zeros(64)
        if (sig_y < 0):# or (sig_y > 18.0):
            return numpy.zeros(64)
        fine_grid = numpy.linspace(0, 8)
        fine_image = numpy.zeros((len(fine_grid), len(fine_grid)))

        a = numpy.cos(theta)**2/(2.0*sig_x**2.0) + numpy.sin(theta)**2./(2.0*sig_y**2.0)
        b = -numpy.sin(2.0*theta)**2/(4.0*sig_x**2.0) + numpy.sin(2.0*theta)**2./(4.0*sig_y**2.0)
        c = numpy.sin(theta)**2/(2.0*sig_x**2.0) + numpy.cos(theta)**2./(2.0*sig_y**2.0)

        indices = range(len(fine_grid))
        
        for x in indices:
            for y in indices:
                fine_image[x][y] = numpy.exp(-(a*(fine_grid[x]-xc)**2.0 +
                    2.0*b*(fine_grid[x]-xc)*(fine_grid[y]-yc) +
                    c*(fine_grid[y]-yc)**2.0))

        coarse_image = numpy.zeros((8,8))
        for x in range(8):
            for y in range(8):
                goodx = numpy.arange(len(fine_grid))[(x <= fine_grid) & (fine_grid < x+1)]
                goody = numpy.arange(len(fine_grid))[(y <= fine_grid) & (fine_grid < y+1)]
                coarse_image[x][y] = numpy.sum(fine_image[goodx][:,goody])

        cutoff = numpy.median(coarse_image) + numpy.std(coarse_image)
        bright = coarse_image > cutoff
        dark = coarse_image <= cutoff
        coarse_image[bright] = 1.0
        coarse_image[dark] /= cutoff
        return coarse_image.ravel()

    def fitElipse(self, cutout):
        """
        cutout should be 8x8 grid of pixels
        """
        sig_x = 2.4
        e = 0.8
        angle = 10.0
        xcenter = 3.5
        ycenter = 3.5

        errfunc = lambda p, y : numpy.abs(self.Elipse(p) - y)
        coeffs = [sig_x, e, angle, xcenter, ycenter]
        pfit = optimize.leastsq(errfunc, coeffs, args=(cutout))

        return pfit

    def extractCutout(self, x, y):
        xcenter = int(numpy.round(self.xcenters[x]))
        ycenter = int(numpy.round(self.ycenters[y]))

        cutout = self.subtracted[ycenter-4:ycenter+4, xcenter-4:xcenter+4]

        return cutout.copy()

    def twirl(self, ax1=None, ax2=None):
        xp = numpy.arange(0, 8, 1.0)
        yp = numpy.arange(0, 8, 1.0)
        grid_x = numpy.linspace(0, 7)
        grid_y = numpy.linspace(0, 7)
        grid = numpy.meshgrid(grid_x, grid_y)
        angles = numpy.arange(0, 360, 45.0)
        for i in range(self.nx):
            for j in range(self.ny):
                x = self.xcenters[i]
                y = self.ycenters[j]
                cutout = self.extractCutout(i, j)
                if not(ax1 == None):
                    ax1.clear()
                    ax2.matshow(cutout)
                    ax1.figure.show()
                if numpy.sum(cutout) > 400:
                    f = scipy.interpolate.interp2d(xp, yp, cutout, kind='cubic')
                    for theta in angles:
                        rotated = rotate(grid, theta)
                        print("%s" % rotated)
                        raw_input()


    def fitElipses(self, ax1=None, ax2=None, ax3=None):
        sig_x = []
        sig_y = []
        angle = []
        amplitude = []
        xc = []
        yc = []
        for i in range(self.nx):
            for j in range(self.ny):
                x = self.xcenters[i]
                y = self.ycenters[j]
                cutout = self.extractCutout(i, j)
                #print x, y, numpy.sum(cutout)
                if not(ax1 == None):
                    ax1.clear()
                    ax1.matshow(cutout, vmin = 0, vmax = 1)
                    ax1.figure.show()
                if numpy.sum(cutout) > 400:
                    cutoff = numpy.median(cutout)+ numpy.std(cutout)
                    good = cutout > cutoff
                    bad = cutout <= cutoff
                    cutout[good] = 1.0
                    cutout[bad] /= cutoff
                    fit, result = self.fitElipse(cutout.ravel())
                    if ax2 != None:
                        ax2.clear()
                        ax2.matshow(self.Elipse(fit).reshape(8,8), vmin = 0, vmax=1)
                        ax2.figure.show()
                        print("%s%" % fit)
                        raw_input()
                    sig_x.append(fit[0])
                    sig_y.append(fit[0]*fit[1])
                    angle.append(fit[2])
                    xc.append(fit[3]-3.5+x)
                    yc.append(fit[3]-3.5+y)
                    print("%.3f, %.3f" %( sig_y[-1], angle[-1]))
        self.sig_x = numpy.array(sig_x)
        self.sig_y = numpy.array(sig_y)
        self.angle = numpy.array(angle)
        self.amplitude = numpy.array(amplitude)
        self.xc = numpy.array(xc)

    def findCentroids(self, ax = None):
        xc = []
        yc = []
        for i in range(self.nx):
            for j in range(self.ny):
                x = self.xcenters[i]
                y = self.ycenters[j]
                cutout = self.extractCutout(i, j)
                if numpy.sum(cutout) > 400:
                    centroid = scipy.ndimage.measurements.center_of_mass(cutout)
                    xc.append(centroid[1] - len(cutout[0])/2.0+x)
                    yc.append(centroid[0] - len(cutout)/2.0+y)
                    if not(ax == None):
                        ax.clear()
                        ax.matshow(cutout)
                        ax.scatter(centroid[1], centroid[0], size=50)
                        ax.figure.show()
                        raw_input()

        self.xc = numpy.array(xc)
        self.yc = numpy.array(yc)

        #self.residual_x = self.xc % 1
        #self.residual_y = self.yc % 1

        #self.residual_x[self.residual_x > 0.5] -= 1.0
        #self.residual_y[self.residual_y > 0.5] -= 1.0


    def fit_line(self, x, y):
        fitfunc = lambda p, x : p[0]+(x*p[1])
        errfunc = lambda p, x, y: numpy.abs(fitfunc(p,x) - y)
        coeffs = [numpy.mean(y), 0.1]
        order = numpy.argsort(x)
        middle_x = numpy.array(x)[order[1:-1]]
        middle_y = numpy.array(y)[order[1:-1]]
        pfit = optimize.leastsq(errfunc, coeffs, args=(middle_x,middle_y) )

        #return numpy.arctan(numpy.abs(pfit[0][1]))*180.0/3.14159262
        return pfit


    def findAngles(self, ax = None):
        xangle = []
        residuals_x = []
        for x in self.xcenters:
            selected = numpy.abs(self.xc - x) < 2.0
            fit = self.fit_line(self.yc[selected], self.xc[selected])
            if not(ax == None):
                ax.plot(fit[0][0]+fit[0][1]*self.yc[selected], 
                        self.yc[selected], color = 'r')
            residuals_x.append(numpy.array([self.yc[selected], 
                fit[0][0]+fit[0][1]*self.yc[selected] - self.xc[selected]]))
            xangle.append(numpy.arctan(fit[0][1])*180.0/3.14159)

        yangle = []
        residuals_y = []
        for y in self.ycenters:
            selected = numpy.abs(self.yc - y) < 2.0
            fit = self.fit_line(self.xc[selected], self.yc[selected])
            if not(ax == None):
                ax.plot(self.xc[selected], fit[0][0]+fit[0][1]*self.xc[selected],
                        color = 'r')
            residuals_y.append(numpy.array([self.xc[selected], 
                fit[0][0]+fit[0][1]*self.xc[selected] - self.yc[selected]]))
            yangle.append(-numpy.arctan(fit[0][1])*180.0/3.14159)

        if not(ax==None):
            ax.scatter(self.xc, self.yc, s=10, c='r')

        self.xangle = xangle
        self.yangle = yangle
        self.residuals_x = residuals_x
        self.residuals_y = residuals_y

    



def twoDgaussian(x, y, center, stdev, A):
    retval = A * (numpy.exp(-(x-center[0])**2.0/stdev[0])*
                  numpy.exp(-(y-center[1])**2.0/stdev[1]))
    print("%.3f %.3f %.3f" % ( center, A, numpy.max(retval), numpy.max(x), numpy.min(x)))
    return retval

class zernikeMode(object):
    """
       Class representing a Zernike mode
    """
    def __init__(self, noll, mag):
        """
        input:  noll - Noll index
                mag - magnitude of Zernike mode - Units?
        """
        i = 0
        n = 0
        m = 0
        while ( i != noll):
            if m == n:
                n += 1
                m = -n
            else:
                m += 2
            i = (n*(n+2) + m)/2.0
        self.n = n
        self.m = m
        self.mag = mag
        """
        if (noll == 1):    # Piston
            self.n = 0
            self.m = 0
        elif (noll == 2):  # Tip
            self.n = 1
            self.m = 1
        elif (noll == 3):  # Tilt
            self.n = 1
            self.m = -1
        elif (noll == 4):  # Defocus
            self.n = 2
            self.m = 0
        elif (noll == 5):
            self.n = 2
            self.m = -2
        elif (noll == 6):
            self.n = 2
            self.m = 2
        elif (noll == 7):
            self.n = 3
            self.m = -1
        elif (noll == 8):
            self.n = 3
            self.m = 1
        elif (noll == 9):
            self.n = 3
            self.m = -3
        elif (noll == 10):
            self.n = 3
            self.m = 3
        elif (noll == 11):
            self.n = 4
            self.m = 0
        elif (noll == 12):
            self.n = 4
            self.m = 2
        elif (noll == 13):
            self.n = 4
            self.m = -2
        elif (noll == 14):
            self.n = 4
            self.m = 4
        elif (noll == 15):
            self.n = 4
            self.m = -4
        else:
            self.n = 0
            self.m = 0
        """

    def zernike_rad(self, rho):
        n=abs(self.n)
        m=abs(self.m)
        if (numpy.mod(n-m, 2) == 1):
            return rho*0.0
        
        wf = rho*0.0
        for k in range((n-m)/2+1):
            wf += rho**(n-2.0*k) * (-1.0)**k * fac(n-k) / ( fac(k) * fac( (n+m)/2.0 - k ) * fac( (n-m)/2.0 - k ) )
        
        return wf

    def setMag(self, mag):
        self.mag = mag

    def zernike(self, rho, phi, norm=False):
        nc = self.mag
        if (norm):
            nc = (2*(self.n+1)/(1+(self.m==0)))**0.5
        if (self.m > 0): return nc*self.zernike_rad(rho) * numpy.cos(self.m*phi)
        if (self.m < 0): return nc*self.zernike_rad(rho) * numpy.sin(-self.m*phi)
        return nc*self.zernike_rad(rho)

class pupil( object ):
    def __init__(self, x, y, innerRadius, outerRadius):
        self.x = x
        self.y = y
        self.inrad2 = innerRadius**2.0
        self.outrad2 = outerRadius**2.0

    def calculateFlux(self, x, y):
        pythag = (self.x-x)**2.0 + (self.y - y)**2.0
        if ( (pythag < self.outrad2) & (pythag > self.inrad2)):
            return 1.0
        else:
            return 0.0

    def calculateApertureIllumination(self, corner, size):
        area = 0.0
        illuminated = 0.0
        for x in numpy.linspace(corner[0], corner[0]+size, 50):
            for y in numpy.linspace(corner[1], corner[1]+size, 50):
                area += 1.0
                pythag = (self.x-x)**2.0 + (self.y - y)**2.0
                if ( (pythag < self.outrad2) & (pythag > self.inrad2)):
                    illuminated += 1.0
        return illuminated/area

    def getDecenter(self):
        return self.x, self.y

    def setDecenter(self, x, y):
        self.x = x
        self.y = y

class deformableMirror( object ):
    def __init__(self, parent=None, IFFile=None):
        self.nActuators = 60
        if IFFile != None:
            self.influenceFunctions = pyfits.getdata(IFFile)
        else:
            self.influenceFunctions = pyfits.getdata(parent.datadir+'IF_cube.fits')
        self.influenceFunctions[numpy.isnan(self.influenceFunctions)] = 0.0
        nx = len(self.influenceFunctions[0][0])
        ny = len(self.influenceFunctions[0])
        xcoords = numpy.linspace(-888.0, 888.0, nx)
        ycoords = numpy.linspace(-888.0, 888.0, ny)
        self.interpFunctions = []
        self.actuatorPositions = numpy.zeros(len(self.influenceFunctions))
        for inflfunc in self.influenceFunctions:
            self.interpFunctions.append(interp.interp2d(xcoords, ycoords,
                            inflfunc, kind='cubic'))

        
    def setMirror(self, actPos):
        self.actuatorPositions = actPos

    def calcPosition(self, x, y):
        retval = 0.0
        for IF, amp in zip(self.interpFunctions, self.actuatorPositions):
            retval += amp*IF(x, y)[0]
        return retval

class derotator( object ):
    def __init__(self, parent):
        self.angle = 0.0
        self.parent = parent

    def setAngle(self, angle):
        self.angle = angle
        self.sine = numpy.sin(self.angle)
        self.cosine = numpy.cos(self.angle)

    def getMirrorPosition(self, x, y):
        newx = x*self.cosine-y*self.sine
        newy = x*self.sine+y*self.cosine
        return newx, newy

class WFS ( object ):
    def __init__(self, wavelength = 1800.0, beamSize = 1776.0, angle = 0.0):
        self.beamSize = beamSize
        self.wavelength = wavelength
        self.centObscScale = 1.116/8.00
        self.datadir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))+'/data/'
        self.detector = detector(self, beamSize = beamSize)
        self.lenslet = lensletArray(self, angle = angle)
        self.wavefront = waveFront(self, nZern=50)
        self.pupil = pupil(0.0, 0.0, innerRadius=self.beamSize/2.0*self.centObscScale,
                outerRadius=self.beamSize/2.0)
        self.DM = deformableMirror(self)
        self.derotator = derotator(self)
        self.centroids = []
        self.WFS_Frame = WFS_Frame()

    def setupInstrument(self, zern=None, pupil=None, actuatorPokes=None,
            derotAngle=0.0, lensletAngle=0.0):
        """
        Generates an image seen by the detector of a wavefront described by
        the zernike coefficients in zern

        zern = [tip, tilt, defocus, astig1, astig2]
        pupil = [x, y]
        actuatorPokes = list of 60 actuator positions
        """
        self.derotator.setAngle(numpy.deg2rad(derotAngle))
        self.lenslet.setClockingAngle(lensletAngle)
        if pupil != None:
            self.pupil.setDecenter(pupil[0], pupil[1])
        self.wavefront.setZern(zern)
        if actuatorPokes != None:
            self.DM.setMirror(actuatorPokes)

    def expose(self):
        self.centroids.append(self.detector.expose(Delta=True))


    def plot(self):
        self.WFS_Frame.gradients = self.centroids[-1]
        self.WFS_Frame.plotGradients()

    def calcWaveFront(self, x, y):
        wave = self.wavefront.calcWaveFront(x, y)
        #rotatedPosition = self.derotator.getMirrorPosition(x, y)
        #mirror = self.DM.calcPosition(rotatedPosition[0], rotatedPosition[1])
        #return wave + mirror
        return wave

class detector( object ):
    """
        The Detector class allows us to simulate the detector in the cryostat.
        Some of the numbers are most likely wrong, but with tweaking we should
        be able to simulate the effects of simple Zernike aberations of the 
        wave front on the spots on the detector.
    """
    def __init__(self, parent, beamSize=1776.0):
        self.parent = parent
        # Calculates the relative size of the central obscuration (M2)
        self.scramblingMap = pyfits.getdata(
                            self.parent.datadir+"scramblemap.fits")
        self.unscramblingMap = pyfits.getdata(
                            self.parent.datadir+"unscramblemap.fits")
        self.windowMap = pyfits.getdata(
                            self.parent.datadir+"windowmap.fits")
        self.SLsubapmap = self.parent.datadir+"LoopDisplaySrv.SUBAP_MAP.fits"
        self.beamSize = beamSize
        self.pixPerSubAp = 8
        self.nx = 72
        self.ny = 72
        self.readoutMode = "8x8"
        self.spacing = 24.0 #microns  Is this value correct?
        self.xpix = (numpy.arange(self.nx)-self.nx/2.0)*self.spacing
        self.ypix = (numpy.arange(self.ny)-self.ny/2.0)*self.spacing
        self.stdev = (8.0*self.spacing, 8.0*self.spacing)
        self.z = []
        self.frames = []
        self.centroids = []
        #self.weightingMap = self.makeWeightingMap()

    def scrambleFrame(self):
        scrambledFrame = numpy.zeros(6912)
        flatframe = self.z[-1].ravel()
        for y, i in zip(flatframe, self.scramblingMap):
            scrambledFrame[i]=y

        self.frames.append(scrambledFrame)
        

    def makeRamp(self):
        z = numpy.zeros((self.ny, self.nx))
        k = 0
        for i in range(self.nx):
            for j in range(self.ny):
                z[i][j] = k
                k += 1
        self.z.append(z)
        self.scrambleFrame()
        self.centroids.append([])

    
    def expose(self):
        debug = False
        centroids = []
        subsamp = 30.0
        FWHM_i = 1.1 * subsamp # FWHM in Pixel space
        FWHM_k = 0.88*self.pixPerSubAp*subsamp/FWHM_i
        location = int(round((self.pixPerSubAp*subsamp - FWHM_k)/2.0))
        delta = self.parent.lenslet.spacing/2.0
        nPixPoints = self.pixPerSubAp*subsamp + 1
        z = numpy.zeros((self.ny, self.nx))
        ptsOnDetector = numpy.linspace(-nPixPoints/2.0*self.spacing, 
                nPixPoints*self.spacing/2.0, nPixPoints)
        gridx, gridy = numpy.meshgrid(ptsOnDetector, ptsOnDetector)
        gridx /= subsamp
        gridy /= subsamp
        totalFlux = 0.0
        if debug:
            fig = pyplot.figure(0)
            fig.clear()
            ax1 = fig.add_axes([0.1, 0.1, 0.4, 0.4])
            ax2 = fig.add_axes([0.1, 0.5, 0.4, 0.4])
            ax3 = fig.add_axes([0.5, 0.1, 0.4, 0.4])
            ax4 = fig.add_axes([0.5, 0.5, 0.4, 0.4])
            extract = []
            count = 1
        for coord in self.parent.lenslet.coordinates:
            flux = numpy.zeros((nPixPoints, nPixPoints), dtype=numpy.complex)
            #intensity = numpy.zeros((nPixPoints, nPixPoints))
            #phase = numpy.zeros((nPixPoints, nPixPoints))
            i = location
            for x in numpy.linspace(coord[0][0]-delta, coord[0][0]+delta, 
                    num=int(round(FWHM_k))):
                j = location
                for y in numpy.linspace(coord[0][1]-delta, coord[0][1]+delta, 
                        num=int(round(FWHM_k))):
                    flux[j][i] = numpy.complex(self.parent.pupil.calculateFlux(x, y),
                            self.parent.calcWaveFront(x, y))
                    #intensity[i][j] = self.pupil.calculateFlux(x, y)
                    #phase[i][j] = self.calcWaveFront(x, y)
                    totalFlux += 1
                    j += 1
                i += 1
            flux[flux.real == 0.0] = numpy.complex(0.0, 0.0)
            sig = abs(flux) > 0
            flux[sig] -= numpy.complex(0, numpy.mean(flux.imag[sig]))
            image = fftpack.fft2(flux)
            image = fftpack.fftshift((image*image.conjugate()).real)
            if debug:
                if count in [8, 10, 15]:
                    extract.append(image)
                    count += 1
            xc = coord[0][0]+gridx
            yc = coord[0][1]+gridy
            inviewx = (self.xpix > numpy.min(xc)) & (self.xpix <= numpy.max(xc))
            inviewy = (self.ypix > numpy.min(yc)) & (self.ypix <= numpy.max(yc))
            weightX = 0.0
            weightY = 0.0
            denom = 0.0
            for i in range(self.nx):
                if inviewx[i]:
                    for j in range(self.ny):
                        if inviewy[j]:
                            fp = scipy.where( (xc >= self.xpix[i]) & (xc < self.xpix[i]+self.spacing+0.5) & (yc >= self.ypix[j]) & (yc < self.ypix[j]+self.spacing+0.5))
                            z[j][i] = numpy.sum(image[fp])
                            weightX += z[j][i]*i
                            weightY += z[j][i]*j
                            denom += z[j][i]
                            #print i, j, z[i][j], coord
                    #raw_input()
            centroids.append([weightX/denom, weightY/denom])

        if debug:
            ax1.matshow(extract[0])
            ax2.matshow(extract[1])
            ax3.matshow(extract[2])
            ax4.matshow(extract[2]-extract[0])
            fig.show()
            print("%.3f" % numpy.max(extract[0]))
            print("%.3f" % numpy.max(extract[1]))
            print("%.3f" % numpy.max(extract[2]))
            print("%.3f" % numpy.max(extract[2]-extract[0]))
            input()
        self.z.append(z)
        self.centroids.append(numpy.array(centroids))
        self.scrambleFrame()

        return self.centroids[-1]


    def measureCentroids(self, frame):
        image = frame.reshape(72,72)
        centroids = []
        for x in range(9):
            for y in range(9):
                if self.parent.lenslet.SLapertureMap[x][y] != 0:
                    ystart = y*8
                    xstart = x*8
                    subaperture = image[ystart:ystart+8, xstart:xstart+8].copy()
                    centroids.append(scipy.ndimage.measurements.center_of_mass(
                        subaperture))
        return numpy.array(centroids)

    def calculateCentroids(self, zern=None, actuatorPokes=None):
        """
            Calcualates the locations of the centroids under the given 
            Zernike coefficients
        """
        if zern != None:
            self.wavefront.setZern(zern)
        if actuatorPokes != None:
            self.DM.setMirror(actuatorPokes)
        dx = 10.0   # Microns
        dy = 10.0   # Microns
        centroids = []
        intensities = []
        DCx, DCy = self.pupil.getDecenter()  # Decenter of Pupil
        for c in self.lenslet.coordinates:
            # Calculates the partial derivatives
            zxp = self.calcWaveFront(c[0]+DCx+dx, c[1]+DCy)
            zxm = self.calcWaveFront(c[0]+DCx-dx, c[1]+DCy)
            zyp = self.calcWaveFront(c[0]+DCx, c[1]+DCy+dy)
            zym = self.calcWaveFront(c[0]+DCx, c[1]+DCy-dy)
            delx = (zxp - zxm)/(2)
            dely = (zyp - zym)/(2)

            # Computes the normal vector to the surface
            normalx = -delx*dy
            normaly = -dely*dx
            normalz = dx*dy

            #Calculates the shift in microns on the detector
            theta_x = scipy.arctan2(normalx, normalz)
            theta_y = scipy.arctan2(normaly, normalz)
            shift_x = scipy.tan(theta_x)*self.lenslet.fl
            shift_y = scipy.tan(theta_y)*self.lenslet.fl
            centroids.append([c[0]+shift_x, c[1]+shift_y])

            intensities.append(self.pupil.calculateApertureIllumination(
                               [c[0]-self.lenslet.spacing/2.0, 
                               c[1]-self.lenslet.spacing/2.0],
                               self.lenslet.spacing))

        return numpy.array(centroids), numpy.array(intensities)

    def makeWeightingMap(self):
        weightingMap = numpy.zeros((72,72),dtype=numpy.float32)

        for coord in self.lenslet.coords:
            x = coord[0]
            y = coord[1]
            print ("Not done yet!")
                    
        
    def calcWaveFront(self, x, y):
        wave = self.wavefront.calcWaveFront(x, y)
        #rotatedPosition = self.derotator.getMirrorPosition(x, y)
        #mirror = self.DM.calcPosition(rotatedPosition[0], rotatedPosition[1])
        #return wave + mirror
        return wave

    def saveRawFrames(self, filename):
        """
        Saves the raw (unscrambled) frames to a data file
        """

        self.z = numpy.array(self.z)
        hdu=pyfits.PrimaryHDU(self.z)
        hdu.writeto(filename, clobber=True)

    def saveFrames(self, filename):
        """
        Saves the frames to a SPARTA-readable data file
        """
        self.frames = numpy.array(self.frames)
        hdu = pyfits.PrimaryHDU(self.frames)
        hdu.scale('int16', bzero=32768, bscale=1)
        hdu.writeto(filename, clobber=True)

    def saveCentroids(self, filename):
        """
        Saves the calculated centroids to a fits file
        """
        self.centroids = numpy.array(self.centroids)
        hdu = pyfits.PrimaryHDU(self.centroids)
        hdu.writeto(filename, clobber=True)

class lensletArray( object ):
    """
    This class simulates the lenslet array
    """
    def __init__(self, parent=None, spacing=192.0, fl=2095.0, angle=0.0):
        """
            Spacing - spacing between adjacent lenslets (in microns)
            fl - focal length of individual lenslet (in microns)
        """
        if parent:
            self.parent = parent
            self.SLapertureMap = pyfits.getdata(self.parent.detector.SLsubapmap)
        
        else:
            self.SLapertureMap = [[False,False,True,True,True,True,True,False,False],
               [False, True, True, True, True, True, True, True, False],
               [True, True, True, True, True, True, True, True, True],
               [True, True, True, True, True, True, True, True, True],
               [True, True, True, True, False, True, True, True, True],
               [True, True, True, True, True, True, True, True, True],
               [True, True, True, True, True, True, True, True, True],
               [False, True, True, True, True, True, True, True, False],
               [False, False, True, True, True, True, True, False, False]]
        self.spacing = spacing   # Lenslet Array Spacing in microns
        self.fl = fl
        self.angle = numpy.deg2rad(angle)
        self.calculateCentroids()

    def setClockingAngle(self, angle):
        self.angle = numpy.deg2rad(angle)
        self.calculateCentroids()

    def calculateCentroids(self):
        coords = []

        for i in range(9):
            for j in range(9):
                if self.SLapertureMap[i][j]:
                    y = (i-4)*self.spacing
                    x = (j-4)*self.spacing
                    coords.append([(x*numpy.cos(self.angle)-y*numpy.sin(self.angle),x*numpy.sin(self.angle)+y*numpy.cos(self.angle))])

        self.coordinates = coords


class waveFront( object ):
    """
    This object describes the wavefront as it hits the lenslet array
    """
    def __init__(self, parent, beamSize = 1776.0, nZern = 12):
        self.parent = parent
        self.beamSize = beamSize
        
        if self.parent == None:
            self.wavelength = 2200.0
        else:
            self.wavelength = self.parent.wavelength
        self.zernikes = []
        self.nZern = nZern
        for i in numpy.arange(2, self.nZern):
            self.zernikes.append(zernikeMode(i, 0.0))

        #self.tip = zernikeMode(2, 0.00)
        #self.tilt = zernikeMode(3, 0.00)
        #self.defocus = zernikeMode(4, 0.00)
        #self.astig1 = zernikeMode(5, 0.00)
        #self.astig2 = zernikeMode(6, 0.00)
    
    def setZern(self, zern):
        """
            Sets the magnitudes of the Zernike components.
        """
        nzern = numpy.min([len(zern), self.nZern])
        for i in range(nzern):
            self.zernikes[i].setMag(zern[i]*2.0*numpy.pi/self.wavelength)
        #for mag, z in zip(zern,
        #        [self.tip, self.tilt, self.defocus, self.astig1, self.astig2]):
        #    z.setMag(mag*2.0*numpy.pi/self.parent.wavelength)

    def calcWaveFront(self, x, y):
        """
        Sums through the different Zernike components at a particular location
        on the wavefront (x, y) to find the local zernike magnitude.
        """
        rho = (x**2.0 + y**2.0)**(0.5)/(self.beamSize/2.0)
        phi = numpy.arctan2(y, x)
        value = 0.0
        for zern in self.zernikes:
            value += zern.zernike(rho, phi)
        return value

    def pokeActuator(self, x, y, InflFunc):
        """
            interpolates the value of an actuator poke at a certain x and y
        """
        return InflFunc(x, y)


