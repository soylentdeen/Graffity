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
from scipy.ndimage.filters import maximum_filter
from scipy.ndimage.morphology import generate_binary_structure, binary_erosion

class PSF( object ):
    def __init__(self, sizeInPix=20, lam=1.6, dlam=0.3, pscale=17.0, M1=8.0, M2=1.3, nLambdaSteps=9):
        self.sizeInPix = sizeInPix
        self.lam = lam*1.0e-6
        self.dlam = dlam*1.0e-6
        self.pscale = pscale * 2.0 *numpy.pi/360.0/60.0**2.0/1.0e3
        self.M1 = M1
        self.M2 = M2
        self.nLambdaSteps = nLambdaSteps
        self.fmax = self.M1 * self.pscale * self.sizeInPix/self.lam

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

    def generateOTF(self):
        total = numpy.zeros([self.sizeInPix, self.sizeInPix])
        for k in numpy.arange(self.nLambdaSteps):
            l = self.lam - self.dlam*(k - self.nLambdaSteps/2.0)/(self.nLambdaSteps - 1)
            fc = self.fmax*self.lam/l

            for i in range(self.sizeInPix):
                for j in range(self.sizeInPix):
                    y = j - self.sizeInPix/2.0 + 0.5
                    x = i - self.sizeInPix/2.0 + 0.5
                    r = numpy.sqrt(x**2.0 + y**2.0)
                    f = r/fc
                    if f < 1:
                        if r < 0.1:
                            total[i, j] += 1.0/nLambdaSteps
                        else:
                            total[i, j] += (self.TelOTF(f, self.M2/self.M1)*self.Sinc(numpy.pi*x/self.sizeInPix)
                                                       *self.Sinc(numpy.pi*y/self.sizeInPix))/self.nLambdaSteps
                    else:
                        total[i, j] += 0.0
        self.OTF = total
        self.PSF = numpy.fft.fftshift(numpy.abs(numpy.fft.fft2(self.OTF)))
        self.PSF = self.normalize(self.PSF)

    def getPSF(self):
        return self.PSF

    def normalize(self, image):
        return image/numpy.max(image)

    def calcStrehl(self, cutout):
        factor = numpy.sum(self.PSF) / numpy.sum(cutout)
        return numpy.max(cutout)*factor
        


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
        


class CircularBuffer( object ):
    def __init__(self, df='', CDMS_BaseDir = '', CDMS_ConfigDir='', S2M = None, ModalBasis = None, Z2DM = None, S2Z = None, HOIM = None, CM=None, TT2HO=None,
                 DM2Z=None, TTM2Z=None, loopRate=500.0, RTC_Delay=0.5e-3):
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
        else:
            self.TT2HO = None
        if CM != None:
            self.CM = pyfits.getdata(CM, ignore_missing_end=True)
            self.CM /= self.controlGain
            self.CM[:60,:] += self.TT2HO.dot(self.CM[60:,:])
        else:
            self.CM = None
        if TTM2Z != None:
            self.TTM2Z = pyfits.getdata(self.CDMS_BaseDir+self.CDMS_ConfigDir+TTM2Z,
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
            print 'OffControl: %E' % self.OffControl
            print 'WFE: %E' % self.WFE
            print 'WFSError: %E' % self.WFSError
            print 'Strehl: %f' % self.Strehl
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
            print ut
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

                #self.x[
                postageStamp.append(numpy.median(numpy.array(PS), axis=0))
                speckleStamp.append(numpy.array(PS))
            self.postageStamps[ut] = numpy.array(postageStamp)
            self.speckleStamps[ut] = numpy.array(speckleStamp)


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
        self.mag = mag
        if (noll == 2):
            self.n = 1
            self.m = 1
        elif (noll == 3):
            self.n = 1
            self.m = -1
        elif (noll == 4):
            self.n = 2
            self.m = 0
        elif (noll == 5):
            self.n = 2
            self.m = -2
        elif (noll == 6):
            self.n = 2
            self.m = 2
        else:
            self.n = 0
            self.m = 0

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
    def __init__(self, parent):
        self.nActuators = 60
        self.influenceFunctions = pyfits.getdata(parent.datadir+'IF_cube.fits')
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
        self.wavefront = waveFront(self)
        self.pupil = pupil(0.0, 0.0, innerRadius=self.beamSize/2.0*self.centObscScale,
                outerRadius=self.beamSize/2.0)
        self.DM = deformableMirror(self)
        self.derotator = derotator(self)
        self.centroids = []

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
        self.centroids.append(self.detector.expose())

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
            for x in numpy.linspace(coord[0]-delta, coord[0]+delta, 
                    num=int(round(FWHM_k))):
                j = location
                for y in numpy.linspace(coord[1]-delta, coord[1]+delta, 
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
            xc = coord[0]+gridx
            yc = coord[1]+gridy
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


    def calculateCentroids(self, zern, actuatorPokes):
        """
            Calcualates the locations of the centroids under the given 
            Zernike coefficients
        """
        self.wavefront.setZern(zern)
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
                    coords.append((x*numpy.cos(self.angle)-y*numpy.sin(self.angle),
                        x*numpy.sin(self.angle)+y*numpy.cos(self.angle)))

        self.coordinates = coords


class waveFront( object ):
    """
    This object describes the wavefront as it hits the lenslet array
    """
    def __init__(self, parent, beamSize = 1776.0, nZern = 12):
        self.parent = parent
        self.beamSize = beamSize
        
        self.zernikes = []
        self.nZern = nZern
        for i in numpy.arange(2, self.nZern):
            self.zernikes.append([i, 0.0])

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
            self.zernikes[i].setMag(zern[i]*2.0*numpy.pi/self.parent.wavelength)
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


