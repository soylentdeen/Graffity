import scipy
import numpy
from matplotlib import pyplot
import Graffity
import glob
import time

fig0 = pyplot.figure(0)
fig0.clear()
ax0 = fig0.add_axes([0.1, 0.1, 0.8, 0.8])

CIAO_dir = '/home/cdeen/Data/CIAO/2016-09-21/CIAO'#/CIAO1/DATA_LOGGER-000038/'
CDMS_Dir = '/home/cdeen/Code/CIAO/SPARTA/SPARTA_CIAO/CONFIG/spcicfg/config/RTCDATA/CIAO/DEFAULT/'
CDMS_BaseDir = '/home/cdeen/Code/CIAO/SPARTA/SPARTA_CIAO/CONFIG/spcicfg'
loopRate = 500 #Hz
RTC_Delay = 0.5e-3

Seeing = {}
Strehl = {}
Tau0 = {}
Time = {}

for ciao in [1,2,3,4]:
    CIAO_datadir = CIAO_dir+str(ciao)+'/'
    CDMSN_Dir = CDMS_BaseDir + str(ciao)+'/config/RTCDATA/CIAO/DEFAULT/'
    seeing = []
    strehl = []
    tau0 = []
    t = []
    for dp, dn, fn in glob.os.walk(CIAO_datadir):
        if len(dn) == 0:
            CB =Graffity.CircularBuffer(dp+'/CIAO_LOOP_0001.fits',S2M=dp+'/RecnOptimiser.S2M_0001.fits', 
                 ModalBasis=CDMSN_Dir+'RecnOptimiser.ModalBasis.fits',Z2DM=CDMSN_Dir+'RecnOptimiser.Z2DM.fits', 
                 CM=dp+'/Recn.REC1.CM_0001.fits', HOIM=dp+'/RecnOptimiser.HO_IM_0001.fits',
                 TT2HO=CDMSN_Dir+'RecnOptimiser.TT2HO.fits', DM2Z=CDMSN_Dir+'RecnOptimiser.DM2Z.fits',
                 TTM2Z=CDMS_Dir+'RecnOptimiser.TTM2Z.fits', loopRate=loopRate, RTC_Delay=RTC_Delay)

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


            print("CIAO #%d"% ciao)
            print("Seeing : %.3f" % (CB.Seeing*CB.Arcsec))
            print("Strehl : %.3f" % CB.Strehl)
            print("Tau0   : %.4f" % CB.Tau0)
            seeing.append(CB.Seeing*CB.Arcsec)
            strehl.append(CB.Strehl)
            tau0.append(CB.Tau0)
            t.append(time.mktime(time.strptime(CB.header.get('ESO TPL START'), '%Y-%m-%dT%H:%M:%S')))
    Seeing[ciao] = numpy.array(seeing)
    Strehl[ciao] = numpy.array(strehl)
    tau0[ciao] = numpy.array(tau0)
    Time[ciao] = numpy.array(t)


ax0.clear()
ax0.scatter(Time[1], Seeing[1], color = 'b')
ax0.scatter(Time[2], Seeing[2], color = 'g')
ax0.scatter(Time[3], Seeing[3], color = 'r')
ax0.scatter(Time[4], Seeing[4], color = 'y')
fig0.show()
raw_input()

ax0.clear()
ax0.scatter(Time[1], Strehl[1], color = 'b')
ax0.scatter(Time[2], Strehl[2], color = 'g')
ax0.scatter(Time[3], Strehl[3], color = 'r')
ax0.scatter(Time[4], Strehl[4], color = 'y')
fig0.show()

raw_input()
ax0.clear()
ax0.scatter(Seeing[1], Strehl[1], color = 'b')
ax0.scatter(Seeing[2], Strehl[2], color = 'g')
ax0.scatter(Seeing[3], Strehl[3], color = 'r')
ax0.scatter(Seeing[4], Strehl[4], color = 'y')
fig0.show()

