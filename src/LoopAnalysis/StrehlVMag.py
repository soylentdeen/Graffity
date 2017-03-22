import scipy
import numpy
from matplotlib import pyplot
import Graffity
import glob
import time

fig0 = pyplot.figure(0)
fig0.clear()
ax0 = fig0.add_axes([0.1, 0.1, 0.8, 0.8])

#CIAO_dir = '/home/cdeen/Data/CIAO/2016-09-21/CIAO'#/CIAO1/DATA_LOGGER-000038/'
#CIAO_dir = '/media/cdeen/My Passport/GRAVITY/CIAO#1/'
CIAO_dir = '/media/cdeen/My Passport/GRAVITY/CIAO#1/Paranal/2016-04-03/Perf_vs_Magn/'
CDMS_Dir = '/data/cdeen/CIAO_Commissioning/spcicfg/config/RTCDATA/CIAO/DEFAULT/'
CDMS_BaseDir = '/data/cdeen/CIAO_Commissioning/spcicfg'
CDMS_ConfigDir = '/config/RTCDATA/CIAO/DEFAULT/'
loopRate = 500 #Hz
RTC_Delay = 0.5e-3

CIAO_datadir = CIAO_dir
seeing = []
strehl = []
rms = []
tau0 = []
t = []
flux = []
#for dp in glob.glob(CIAO_datadir+'DATA_LOGGER*'):
#    if (dp.find('DATA_LOGGER') != -1):
for dp, dn, fn in glob.os.walk(CIAO_datadir):
    if (len(dn) == 0) & (dp.find('AO_TF') != -1):
        try:
            CB =Graffity.CircularBuffer(dp+'/CIAO_LOOP_0002.fits', CDMS_BaseDir=CDMS_BaseDir, 
                 CDMS_ConfigDir = CDMS_ConfigDir, S2M=dp+'/RecnOptimiser.S2M_0002.fits', 
                 ModalBasis='RecnOptimiser.ModalBasis.fits',Z2DM='RecnOptimiser.Z2DM.fits', 
                 CM=dp+'/Recn.REC1.CM_0002.fits', HOIM=dp+'/RecnOptimiser.HO_IM_0002.fits',
                 TT2HO='RecnOptimiser.TT2HO.fits', DM2Z='RecnOptimiser.DM2Z.fits',
                 TTM2Z='RecnOptimiser.TTM2Z.fits', loopRate=loopRate, RTC_Delay=RTC_Delay)

        
            tm = time.mktime(time.strptime(CB.header.get('ESO TPL START'), '%Y-%m-%dT%H:%M:%S'))/3600.0
            CB.synchronizeData()
            CB.zernikeSpace()
            CB.computePSD(source = 'ZSlopes')
    
            CB.computePSD(source = 'ZCommands')
    
            CB.AORejectionFunction()
            CB.combinePSDs()
            CB.computeKolmogorovCovar()
            CB.zernikePropagation()
            CB.noiseEvaluation()
            CB.seeingEstimation()
            CB.computeSpectralSlope((3,15))
            CB.estimateStrehlRatio()


            print("%.3f" % CB.header.get('ESO AOS LOOP RATE'))
            print("%d" % CB.CIAO_ID)
            print("file: %s " % dp)
            print("Seeing : %.3f" % (CB.Seeing*CB.Arcsec))
            print("Strehl : %.3f" % CB.Strehl)
            print("Tau0   : %.4f" % CB.Tau0)
            seeing.append(CB.Seeing*CB.Arcsec)
            strehl.append(CB.Strehl)
            rms.append(CB.rms)
            tau0.append(CB.Tau0)
            t.append(tm)
            flux.append(numpy.median(CB.Intensities))
        except:
            print asdf
            pass


seeing = numpy.array(seeing)
strehl = numpy.array(strehl)
flux = numpy.array(flux)

#ax0.clear()
#ax0.scatter(Time[1], Seeing[1], color = 'b')
#ax0.scatter(Time[2], Seeing[2], color = 'g')
#ax0.scatter(Time[3], Seeing[3], color = 'r')
#ax0.scatter(Time[4], Seeing[4], color = 'y')
#fig0.show()
#raw_input()

#"""
ax0.clear()
#ax0.scatter(Seeing[1], RMS[1], color = 'b')
#ax0.scatter(Seeing[2], RMS[2], color = 'g')
#ax0.scatter(Seeing[3], RMS[3], color = 'r')
#ax0.scatter(Seeing[4], RMS[4], color = 'y')
ax0.scatter(flux[seeing > 0.2], strehl[seeing > 0.2], color = 'b', label='UT1')
ax0.set_xlabel("Estimated Seeing (\")")
ax0.set_ylabel("Estimated Strehl Ratio")
ax0.legend(loc=1, scatterpoints=1)
pyplot.rcParams["font.size"] = 15.0
fig0.show()
fig0.savefig("StrehlvMag.png")
#"""
