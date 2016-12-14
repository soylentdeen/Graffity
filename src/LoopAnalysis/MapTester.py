import scipy
import numpy
from matplotlib import pyplot
import Graffity
import glob
import time
import astropy.io.fits as pyfits

fig0 = pyplot.figure(0)
fig0.clear()
ax1 = fig0.add_axes([0.1, 0.1, 0.4, 0.4])
ax2 = fig0.add_axes([0.1, 0.5, 0.4, 0.4])
ax3 = fig0.add_axes([0.5, 0.1, 0.4, 0.4])
ax4 = fig0.add_axes([0.5, 0.5, 0.4, 0.4])

CIAO_dir = '/home/cdeen/Data/CIAO/2016-09-21/CIAO'
CDMS_Dir = '/home/cdeen/Code/CIAO/SPARTA/SPARTA_CIAO/CONFIG/spcicfg/config/RTCDATA/CIAO/DEFAULT/'
CDMS_BaseDir = '/home/cdeen/Code/CIAO/SPARTA/SPARTA_CIAO/CONFIG/spcicfg'

#RecnOptimiser.ModalBasis
modalBasis_1 = pyfits.getdata(CDMS_BaseDir+'1/config/RTCDATA/CIAO/DEFAULT/RecnOptimiser.ModalBasis.fits', ignore_missing_end=True)
modalBasis_2 = pyfits.getdata(CDMS_BaseDir+'2/config/RTCDATA/CIAO/DEFAULT/RecnOptimiser.ModalBasis.fits', ignore_missing_end=True)
modalBasis_3 = pyfits.getdata(CDMS_BaseDir+'3/config/RTCDATA/CIAO/DEFAULT/RecnOptimiser.ModalBasis.fits', ignore_missing_end=True)
modalBasis_4 = pyfits.getdata(CDMS_BaseDir+'4/config/RTCDATA/CIAO/DEFAULT/RecnOptimiser.ModalBasis.fits', ignore_missing_end=True)

ax1.matshow(modalBasis_1)
ax2.matshow(modalBasis_2)
ax3.matshow(modalBasis_3)
ax4.matshow(modalBasis_4)
fig0.show()
raw_input()
#RecnOptimiser.Z2DM
ax1.clear()
ax2.clear()
ax3.clear()
ax4.clear()
Z2DM_1 = pyfits.getdata(CDMS_BaseDir+'1/config/RTCDATA/CIAO/DEFAULT/RecnOptimiser.Z2DM.fits', ignore_missing_end=True)
Z2DM_2 = pyfits.getdata(CDMS_BaseDir+'2/config/RTCDATA/CIAO/DEFAULT/RecnOptimiser.Z2DM.fits', ignore_missing_end=True)
Z2DM_3 = pyfits.getdata(CDMS_BaseDir+'3/config/RTCDATA/CIAO/DEFAULT/RecnOptimiser.Z2DM.fits', ignore_missing_end=True)
Z2DM_4 = pyfits.getdata(CDMS_BaseDir+'4/config/RTCDATA/CIAO/DEFAULT/RecnOptimiser.Z2DM.fits', ignore_missing_end=True)

ax1.matshow(Z2DM_1)
ax2.matshow(Z2DM_2)
ax3.matshow(Z2DM_3)
ax4.matshow(Z2DM_4)
fig0.show()
raw_input()

#RecnOptimiser.TT2HO
ax1.clear()
ax2.clear()
ax3.clear()
ax4.clear()
TT2HO_1 = pyfits.getdata(CDMS_BaseDir+'1/config/RTCDATA/CIAO/DEFAULT/RecnOptimiser.TT2HO.fits', ignore_missing_end=True)
TT2HO_2 = pyfits.getdata(CDMS_BaseDir+'2/config/RTCDATA/CIAO/DEFAULT/RecnOptimiser.TT2HO.fits', ignore_missing_end=True)
TT2HO_3 = pyfits.getdata(CDMS_BaseDir+'3/config/RTCDATA/CIAO/DEFAULT/RecnOptimiser.TT2HO.fits', ignore_missing_end=True)
TT2HO_4 = pyfits.getdata(CDMS_BaseDir+'4/config/RTCDATA/CIAO/DEFAULT/RecnOptimiser.TT2HO.fits', ignore_missing_end=True)

ax1.matshow(TT2HO_1)
ax2.matshow(TT2HO_2)
ax3.matshow(TT2HO_3)
ax4.matshow(TT2HO_4)
fig0.show()
raw_input()


#RecnOptimiser.DM2Z
ax1.clear()
ax2.clear()
ax3.clear()
ax4.clear()
DM2Z_1 = pyfits.getdata(CDMS_BaseDir+'1/config/RTCDATA/CIAO/DEFAULT/RecnOptimiser.DM2Z.fits', ignore_missing_end=True)
DM2Z_2 = pyfits.getdata(CDMS_BaseDir+'2/config/RTCDATA/CIAO/DEFAULT/RecnOptimiser.DM2Z.fits', ignore_missing_end=True)
DM2Z_3 = pyfits.getdata(CDMS_BaseDir+'3/config/RTCDATA/CIAO/DEFAULT/RecnOptimiser.DM2Z.fits', ignore_missing_end=True)
DM2Z_4 = pyfits.getdata(CDMS_BaseDir+'4/config/RTCDATA/CIAO/DEFAULT/RecnOptimiser.DM2Z.fits', ignore_missing_end=True)

ax1.matshow(DM2Z_1)
ax2.matshow(DM2Z_2)
ax3.matshow(DM2Z_3)
ax4.matshow(DM2Z_4)
fig0.show()
raw_input()


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
            t.append((time.mktime(time.strptime(CB.header.get('ESO TPL START'), '%Y-%m-%dT%H:%M:%S'))-startTime)/3600.0)
    Seeing[ciao] = numpy.array(seeing)
    Strehl[ciao] = numpy.array(strehl)
    tau0[ciao] = numpy.array(tau0)
    Time[ciao] = numpy.array(t)

#ax0.clear()
#ax0.scatter(Time[1], Seeing[1], color = 'b')
#ax0.scatter(Time[2], Seeing[2], color = 'g')
#ax0.scatter(Time[3], Seeing[3], color = 'r')
#ax0.scatter(Time[4], Seeing[4], color = 'y')
#fig0.show()
#raw_input()

binSize = numpy.mean(numpy.diff(psf_Time))

sr_ut1 = []
dsr_ut1 = []
sr_ut2 = []
dsr_ut2 = []
sr_ut3 = []
dsr_ut3 = []
sr_ut4 = []
dsr_ut4 = []
for i in range(len(psf_Time)):
    if i == 0:
        start = psf_Time[i] - binSize
    else:
        start = (psf_Time[i-1] + psf_Time[i])/2.0
    if i == len(psf_Time)-1:
        stop = psf_Time[i] + binSize
    else:
        stop = (psf_Time[i] + psf_Time[i+1])/2.0

    measurements = (start < Time[1]) & (Time[1] < stop)
    sr_ut1.append(numpy.mean(Strehl[1][measurements]))
    dsr_ut1.append(numpy.std(Strehl[1][measurements]))
    measurements = (start < Time[2]) & (Time[2] < stop)
    sr_ut2.append(numpy.mean(Strehl[2][measurements]))
    dsr_ut2.append(numpy.std(Strehl[2][measurements]))
    measurements = (start < Time[3]) & (Time[3] < stop)
    sr_ut3.append(numpy.mean(Strehl[3][measurements]))
    dsr_ut3.append(numpy.std(Strehl[3][measurements]))
    measurements = (start < Time[4]) & (Time[4] < stop)
    sr_ut4.append(numpy.mean(Strehl[4][measurements]))
    dsr_ut4.append(numpy.std(Strehl[4][measurements]))

"""
#ax0.clear()
ax0.scatter(Time[1], Strehl[1], color = 'b')
ax0.errorbar(psf_Time, sr_ut1, dsr_ut1, color = 'b')
ax0.scatter(Time[2], Strehl[2], color = 'g')
ax0.errorbar(psf_Time, sr_ut2, dsr_ut2, color = 'g')
ax0.scatter(Time[3], Strehl[3], color = 'r')
ax0.errorbar(psf_Time, sr_ut3, dsr_ut3, color = 'r')
ax0.scatter(Time[4], Strehl[4], color = 'y')
ax0.errorbar(psf_Time, sr_ut4, dsr_ut4, color = 'y')

ax0.set_xbound(0, 4.0)
ax0.set_ybound(0.0, 0.9)
ax0.set_xlabel("Hours since start of Observing")
ax0.set_ylabel("Computed Strehl Ratio")
pyplot.rcParams["font.size"] = 20.0
ax0.legend(loc = 4)

fig0.show()
"""
#raw_input()
ax0.clear()
ax0.scatter(Seeing[1], Strehl[1], color = 'b')
ax0.scatter(Seeing[2], Strehl[2], color = 'g')
ax0.scatter(Seeing[3], Strehl[3], color = 'r')
ax0.scatter(Seeing[4], Strehl[4], color = 'y')
fig0.show()

