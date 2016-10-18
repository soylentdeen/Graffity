import scipy
import numpy
from matplotlib import pyplot
import Graffity

fig0 = pyplot.figure(0)
fig0.clear()
ax0 = fig0.add_axes([0.1, 0.1, 0.8, 0.8])

CIAO_datadir = '/data/cdeen/Data/CIAO/Paranal/2016-09-23_4/DATA_LOGGER-003612/'
CDMS_Dir = '/data/cdeen/CIAO_Commissioning/spcicfg/config/RTCDATA/CIAO/DEFAULT/'
CDMS4_Dir = '/data/cdeen/CIAO_Commissioning/spcicfg4/config/RTCDATA/CIAO/DEFAULT/'
loopRate = 500 #Hz
RTC_Delay = 0.5e-3

CB = Graffity.CircularBuffer(CIAO_datadir+'CIAO_LOOP_0001.fits', S2M=CIAO_datadir+'RecnOptimiser.S2M_0001.fits', 
                             ModalBasis=CDMS4_Dir+'RecnOptimiser.ModalBasis.fits', Z2DM=CDMS4_Dir+'RecnOptimiser.Z2DM.fits', 
                             CM=CIAO_datadir+'Recn.REC1.CM_0001.fits', HOIM=CIAO_datadir+'RecnOptimiser.HO_IM_0001.fits',
                             TT2HO=CDMS4_Dir+'RecnOptimiser.TT2HO.fits', DM2Z=CDMS4_Dir+'RecnOptimiser.DM2Z.fits',
                             TTM2Z=CDMS_Dir+'RecnOptimiser.TTM2Z.fits', loopRate=loopRate, RTC_Delay=RTC_Delay)

CB.synchronizeData()
CB.zernikeSpace()
CB.computePSD(source = 'ZSlopes')

"""
ax0.plot(CB.ZPowerFrequencies, CB.ZPowerSlopes)
ax0.set_xscale('log')
ax0.set_yscale('log')
fig0.show()
raw_input()
"""
CB.computePSD(source = 'ZCommands')

"""
ax0.clear()
ax0.plot(CB.ZPowerFrequencies, CB.ZPowerCommands)
ax0.set_xscale('log')
ax0.set_yscale('log')
fig0.show()
raw_input()
"""

CB.AORejectionFunction()
CB.combinePSDs()
CB.computeKolmogorovCovar()
CB.zernikePropagation()
CB.noiseEvaluation()
CB.seeingEstimation()