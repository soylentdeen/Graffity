import scipy
import numpy
from matplotlib import pyplot
import Graffity


CIAO_datadir = '/media/cdeen/My Passport/GRAVITY/CIAO#4/2016-09-20_4/DATA_LOGGER-235016/'
CDMS_Dir = '/home/cdeen/Code/CIAO/SPARTA/SPARTA_CIAO/CONFIG/spcicfg/config/RTCDATA/CIAO/DEFAULT/'
CDMS4_Dir = '/home/cdeen/Code/CIAO/SPARTA/SPARTA_CIAO/CONFIG/spcicfg4/config/RTCDATA/CIAO/DEFAULT/'
loopRate = 500 #Hz
RTC_Delay = 0.5e-3

CB = Graffity.CircularBuffer(CIAO_datadir+'CIAO_LOOP_0001.fits', S2M=CIAO_datadir+'RecnOptimiser.S2M_0001.fits', 
                             ModalBasis=CDMS4_Dir+'RecnOptimiser.ModalBasis.fits', Z2DM=CDMS4_Dir+'RecnOptimiser.Z2DM.fits', 
                             CM=CIAO_datadir+'Recn.REC1.CM_0001.fits', HOIM=CIAO_datadir+'RecnOptimiser.HO_IM_0001.fits',
                             TT2HO=CDMS4_Dir+'RecnOptimiser.TT2HO.fits', DM2Z=CDMS4_Dir+'RecnOptimiser.DM2Z.fits',
                             TTM2Z=CDMS_Dir+'RecnOptimiser.TTM2Z.fits', loopRate=loopRate, RTC_Delay=RTC_Delay)

CB.synchronizeData()
CB.zernikeSpace()