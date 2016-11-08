import scipy
import numpy
from matplotlib import pyplot
import Graffity
import glob
import os
import time

fig0 = pyplot.figure(0)
fig0.clear()
ax0 = fig0.add_axes([0.1, 0.1, 0.8, 0.8])

CDMS_Dir = '/data/cdeen/CIAO_Commissioning/spcicfg/config/RTCDATA/CIAO/DEFAULT/'
CDMS4_Dir = '/data/cdeen/CIAO_Commissioning/spcicfg4/config/RTCDATA/CIAO/DEFAULT/'
loopRate = 500 #Hz
RTC_Delay = 0.5e-3

dataloggers = glob.glob('/data/cdeen/Data/CIAO/Paranal/2016-09-23_4/DATA_LOGGER-*')

Time = []
Seeing = []
Strehl = []
Tau0 = []

for d in dataloggers:
    CIAO_datadir = d+'/' #'/data/cdeen/Data/CIAO/Paranal/2016-09-23_4/DATA_LOGGER-003612/'
    CB = Graffity.CircularBuffer(CIAO_datadir+'CIAO_LOOP_0001.fits', S2M=CIAO_datadir+'RecnOptimiser.S2M_0001.fits', 
                             ModalBasis=CDMS4_Dir+'RecnOptimiser.ModalBasis.fits', Z2DM=CDMS4_Dir+'RecnOptimiser.Z2DM.fits', 
                             CM=CIAO_datadir+'Recn.REC1.CM_0001.fits', HOIM=CIAO_datadir+'RecnOptimiser.HO_IM_0001.fits',
                             TT2HO=CDMS4_Dir+'RecnOptimiser.TT2HO.fits', DM2Z=CDMS4_Dir+'RecnOptimiser.DM2Z.fits',
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

    print("Seeing : %.3f" % (CB.Seeing*CB.Arcsec))
    print("Strehl : %.3f" % CB.Strehl)
    print("Tau0   : %.4f" % CB.Tau0)

    Seeing.append(CB.Seeing*CB.Arcsec)
    Strehl.append(CB.Strehl)
    Tau0.append(CB.Tau0)
    Time.append(time.mktime(time.strptime(CB.header.get('ESO TPL START'), '%Y-%m-%dT%H:%M:%S')))
    del(CB)

Seeing = numpy.array(Seeing)
Strehl = numpy.array(Strehl)
Tau0 = numpy.array(Tau0)
Time = numpy.array(Time)

ax0.clear()
#ax0.scatter(Time, Seeing, color = 'b')
#ax0.scatter(Time, Strehl, color = 'r')
ax0.scatter(Seeing, Strehl, color = 'g')
ax0.set_xlabel("Seeing")
ax0.set_ylabel("Strehl Ratio")
ax0.set_title("Seeing Vs. Strehl - 23 Sept, UT4")
fig0.show()
fig0.savefig("23Sept_UT4.png")
