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
CDMS_ConfigDir = '/config/RTCDATA/CIAO/DEFAULT/'
loopRate = 500 #Hz
RTC_Delay = 0.5e-3

Seeing = {}
Strehl = {}
RMS = {}
Tau0 = {}
Time = {}

startTime = time.mktime(time.strptime('2016-09-21T23:00:00', '%Y-%m-%dT%H:%M:%S'))

PSFStrehl = open('strehlRatios.dat', 'r')
psf_Time = []
psf_ut1 = []
psf_ut2 = []
psf_ut3 = []
psf_ut4 = []
psf_dut1 = []
psf_dut2 = []
psf_dut3 = []
psf_dut4 = []
for line in PSFStrehl.readlines():
    l = line.split()
    psf_Time.append((numpy.int(l[0])- startTime)/3600.0)
    psf_ut1.append(numpy.exp(numpy.log(numpy.float(l[1]))*(1.6/2.2)**2.0))
    psf_dut1.append(numpy.exp(numpy.log(numpy.float(l[2]))*(1.6/2.2)**2.0))
    psf_ut2.append(numpy.exp(numpy.log(numpy.float(l[3]))*(1.6/2.2)**2.0))
    psf_dut2.append(numpy.exp(numpy.log(numpy.float(l[4]))*(1.6/2.2)**2.0))
    psf_ut3.append(numpy.exp(numpy.log(numpy.float(l[5]))*(1.6/2.2)**2.0))
    psf_dut3.append(numpy.exp(numpy.log(numpy.float(l[6]))*(1.6/2.2)**2.0))
    psf_ut4.append(numpy.exp(numpy.log(numpy.float(l[7]))*(1.6/2.2)**2.0))
    psf_dut4.append(numpy.exp(numpy.log(numpy.float(l[8]))*(1.6/2.2)**2.0))

psf_Time = numpy.array(psf_Time)
psf_ut1 = numpy.array(psf_ut1)
psf_dut1 = numpy.array(psf_dut1)
psf_ut2 = numpy.array(psf_ut2)
psf_dut2 = numpy.array(psf_dut2)
psf_ut3 = numpy.array(psf_ut3)
psf_dut3 = numpy.array(psf_dut3)
psf_ut4 = numpy.array(psf_ut4)
psf_dut4 = numpy.array(psf_dut4)


for ciao in [1,2,3,4]:
    CIAO_datadir = CIAO_dir+str(ciao)+'/'
    seeing = []
    strehl = []
    rms = []
    tau0 = []
    t = []
    for dp, dn, fn in glob.os.walk(CIAO_datadir):
        if len(dn) == 0:
            CB =Graffity.CircularBuffer(dp+'/CIAO_LOOP_0001.fits', CDMS_BaseDir=CDMS_BaseDir, 
                 CDMS_ConfigDir = CDMS_ConfigDir, S2M=dp+'/RecnOptimiser.S2M_0001.fits', 
                 ModalBasis='RecnOptimiser.ModalBasis.fits',Z2DM='RecnOptimiser.Z2DM.fits', 
                 CM=dp+'/Recn.REC1.CM_0001.fits', HOIM=dp+'/RecnOptimiser.HO_IM_0001.fits',
                 TT2HO='RecnOptimiser.TT2HO.fits', DM2Z='RecnOptimiser.DM2Z.fits',
                 TTM2Z='RecnOptimiser.TTM2Z.fits', loopRate=loopRate, RTC_Delay=RTC_Delay)

            
            tm = (time.mktime(time.strptime(CB.header.get('ESO TPL START'), '%Y-%m-%dT%H:%M:%S'))-startTime)/3600.0
            if (( 0 < tm ) & ( tm < 24)):
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
    Seeing[ciao] = numpy.array(seeing)
    Strehl[ciao] = numpy.array(strehl)
    RMS[ciao] = numpy.array(rms)
    Tau0[ciao] = numpy.array(tau0)
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

seeing_ut1 = []
dseeing_ut1 = []
seeing_ut2 = []
dseeing_ut2 = []
seeing_ut3 = []
dseeing_ut3 = []
seeing_ut4 = []
dseeing_ut4 = []

slopes_ut1 = []
dslopes_ut1 = []
slopes_ut2 = []
dslopes_ut2 = []
slopes_ut3 = []
dslopes_ut3 = []
slopes_ut4 = []
dslopes_ut4 = []

tau_ut1 = []
dtau_ut1 = []
tau_ut2 = []
dtau_ut2 = []
tau_ut3 = []
dtau_ut3 = []
tau_ut4 = []
dtau_ut4 = []

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
    seeing_ut1.append(numpy.mean(Seeing[1][measurements]))
    dseeing_ut1.append(numpy.std(Seeing[1][measurements]))
    slopes_ut1.append(numpy.mean(RMS[1][measurements]))
    dslopes_ut1.append(numpy.std(RMS[1][measurements]))
    tau_ut1.append(numpy.mean(Tau0[1][measurements]))
    dtau_ut1.append(numpy.std(Tau0[1][measurements]))

    measurements = (start < Time[2]) & (Time[2] < stop)
    sr_ut2.append(numpy.mean(Strehl[2][measurements]))
    dsr_ut2.append(numpy.std(Strehl[2][measurements]))
    seeing_ut2.append(numpy.mean(Seeing[2][measurements]))
    dseeing_ut2.append(numpy.std(Seeing[2][measurements]))
    slopes_ut2.append(numpy.mean(RMS[2][measurements]))
    dslopes_ut2.append(numpy.std(RMS[2][measurements]))
    tau_ut2.append(numpy.mean(Tau0[2][measurements]))
    dtau_ut2.append(numpy.std(Tau0[2][measurements]))

    measurements = (start < Time[3]) & (Time[3] < stop)
    sr_ut3.append(numpy.mean(Strehl[3][measurements]))
    dsr_ut3.append(numpy.std(Strehl[3][measurements]))
    seeing_ut3.append(numpy.mean(Seeing[3][measurements]))
    dseeing_ut3.append(numpy.std(Seeing[3][measurements]))
    slopes_ut3.append(numpy.mean(RMS[3][measurements]))
    dslopes_ut3.append(numpy.std(RMS[3][measurements]))
    tau_ut3.append(numpy.mean(Tau0[3][measurements]))
    dtau_ut3.append(numpy.std(Tau0[3][measurements]))

    measurements = (start < Time[4]) & (Time[4] < stop)
    sr_ut4.append(numpy.mean(Strehl[4][measurements]))
    dsr_ut4.append(numpy.std(Strehl[4][measurements]))
    seeing_ut4.append(numpy.mean(Seeing[4][measurements]))
    dseeing_ut4.append(numpy.std(Seeing[4][measurements]))
    slopes_ut4.append(numpy.mean(RMS[4][measurements]))
    dslopes_ut4.append(numpy.std(RMS[4][measurements]))
    tau_ut4.append(numpy.mean(Tau0[4][measurements]))
    dtau_ut4.append(numpy.std(Tau0[4][measurements]))



sr_ut1 = numpy.array(sr_ut1)
dsr_ut1 = numpy.array(dsr_ut1)
seeing_ut1 = numpy.array(seeing_ut1)
dseeing_ut1 = numpy.array(dseeing_ut1)
slopes_ut1 = numpy.array(slopes_ut1)
dslopes_ut1 = numpy.array(dslopes_ut1)
tau_ut1 = numpy.array(tau_ut1)
dtau_ut1 = numpy.array(dtau_ut1)

sr_ut2 = numpy.array(sr_ut2)
dsr_ut2 = numpy.array(dsr_ut2)
seeing_ut2 = numpy.array(seeing_ut2)
dseeing_ut2 = numpy.array(dseeing_ut2)
slopes_ut2 = numpy.array(slopes_ut2)
dslopes_ut2 = numpy.array(dslopes_ut2)
tau_ut2 = numpy.array(tau_ut2)
dtau_ut2 = numpy.array(dtau_ut2)

sr_ut3 = numpy.array(sr_ut3)
dsr_ut3 = numpy.array(dsr_ut3)
seeing_ut3 = numpy.array(seeing_ut3)
dseeing_ut3 = numpy.array(dseeing_ut3)
slopes_ut3 = numpy.array(slopes_ut3)
dslopes_ut3 = numpy.array(dslopes_ut3)
tau_ut3 = numpy.array(tau_ut3)
dtau_ut3 = numpy.array(dtau_ut3)

sr_ut4 = numpy.array(sr_ut4)
dsr_ut4 = numpy.array(dsr_ut4)
seeing_ut4 = numpy.array(seeing_ut4)
dseeing_ut4 = numpy.array(dseeing_ut4)
slopes_ut4 = numpy.array(slopes_ut4)
dslopes_ut4 = numpy.array(dslopes_ut4)
tau_ut4 = numpy.array(tau_ut4)
dtau_ut4 = numpy.array(dtau_ut4)

"""
#ax0.clear()

A1 = []
A2 = []
A3 = []
A4 = []
for SR1, SR2, SR3, SR4, seeing1, seeing2, seeing3, seeing4, tau1, tau2, tau3, tau4 in zip(sr_ut1, sr_ut2, sr_ut3, sr_ut4, seeing_ut1, seeing_ut2, seeing_ut3, seeing_ut4, tau_ut1, tau_ut2, tau_ut3, tau_ut4):
    A1.append([1.0, SR1, seeing1, tau1])
    A2.append([1.0, SR2, seeing2, tau2])
    A3.append([1.0, SR3, seeing3, tau3])
    A4.append([1.0, SR4, seeing4, tau4])

A1 = numpy.array(A1)
B1 = numpy.linalg.pinv(A1)
fit1 = B1.dot(psf_ut1)

A2 = numpy.array(A2)
B2 = numpy.linalg.pinv(A2)
fit2 = B2.dot(psf_ut2)

A3 = numpy.array(A3)
B3 = numpy.linalg.pinv(A3)
fit3 = B3.dot(psf_ut3)

A4 = numpy.array(A4)
B4 = numpy.linalg.pinv(A4)
fit4 = B4.dot(psf_ut4)

print fit1
print fit2
print fit3
print fit4
ax0.clear()
ax0.scatter(psf_ut1, sr_ut1*fit1[1] + fit1[0] + seeing_ut1*fit1[2] + tau_ut1*fit1[3], color = 'b')
ax0.scatter(psf_ut2, sr_ut2*fit2[1] + fit2[0] + seeing_ut2*fit2[2] + tau_ut2*fit2[3], color = 'r')
ax0.scatter(psf_ut3, sr_ut3*fit3[1] + fit3[0] + seeing_ut3*fit3[2] + tau_ut3*fit3[3], color = 'g')
ax0.scatter(psf_ut4, sr_ut4*fit4[1] + fit4[0] + seeing_ut4*fit4[2] + tau_ut4*fit4[3], color = 'm')
ax0.plot([0.35, 0.55], [0.35, 0.55])
fig0.show()

raw_input()
ax0.clear()

ax0.errorbar(psf_Time, psf_ut1, psf_dut1, color = 'b', label = 'UT1')
ax0.errorbar(psf_Time, psf_ut2, psf_dut2, color = 'g', label = 'UT2')
ax0.errorbar(psf_Time, psf_ut3, psf_dut3, color = 'r', label = 'UT3')
ax0.errorbar(psf_Time, psf_ut4, psf_dut4, color = 'y', label = 'UT4')

ax0.scatter(Time[1], Strehl[1], color = 'b')
ax0.errorbar(psf_Time, sr_ut1, dsr_ut1, color = 'b')
ax0.scatter(psf_Time, sr_ut1*fit1[1]+fit1[0]+seeing_ut1*fit1[2]+tau_ut1*fit1[3], color = 'b')
ax0.scatter(Time[2], Strehl[2], color = 'g')
ax0.errorbar(psf_Time, sr_ut2, dsr_ut2, color = 'g')
ax0.scatter(psf_Time, sr_ut2*fit2[1]+fit2[0]+seeing_ut2*fit2[2]+tau_ut2*fit2[3], color = 'g')
ax0.scatter(Time[3], Strehl[3], color = 'r')
ax0.errorbar(psf_Time, sr_ut3, dsr_ut3, color = 'r')
ax0.scatter(psf_Time, sr_ut3*fit3[1]+fit3[0]+seeing_ut3*fit3[2]+tau_ut3*fit3[3], color = 'r')
ax0.scatter(Time[4], Strehl[4], color = 'y')
ax0.errorbar(psf_Time, sr_ut4, dsr_ut4, color = 'y')
ax0.scatter(psf_Time, sr_ut4*fit4[1]+fit4[0]+seeing_ut4*fit4[2]+tau_ut4*fit4[3], color = 'y')

ax0.set_xbound(0, 4.0)
ax0.set_ybound(0.0, 0.9)
ax0.set_xlabel("Hours since start of Observing")
ax0.set_ylabel("Computed Strehl Ratio")
pyplot.rcParams["font.size"] = 20.0
ax0.legend(loc = 4)

fig0.show()
#"""
#"""
ax0.clear()
#ax0.scatter(Seeing[1], RMS[1], color = 'b')
#ax0.scatter(Seeing[2], RMS[2], color = 'g')
#ax0.scatter(Seeing[3], RMS[3], color = 'r')
#ax0.scatter(Seeing[4], RMS[4], color = 'y')
ax0.scatter(Seeing[1][Seeing[1] > 0.2], Strehl[1][Seeing[1] > 0.2], color = 'b', label='UT1')
ax0.scatter(Seeing[2][Seeing[2] > 0.2], Strehl[2][Seeing[2] > 0.2], color = 'g', label='UT2')
ax0.scatter(Seeing[3][Seeing[3] > 0.2], Strehl[3][Seeing[3] > 0.2], color = 'r', label='UT3')
ax0.scatter(Seeing[4][Seeing[4] > 0.2], Strehl[4][Seeing[4] > 0.2], color = 'y', label='UT4')
ax0.set_xlabel("Estimated Seeing (\")")
ax0.set_ylabel("Estimated Strehl Ratio")
ax0.legend(loc=1, scatterpoints=1)
pyplot.rcParams["font.size"] = 15.0
fig0.show()
fig0.savefig("StrehlvSeeing.png")
#"""
