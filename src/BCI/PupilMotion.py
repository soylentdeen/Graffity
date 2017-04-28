import Graffity
import CIAO_DatabaseTools
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D
import numpy
import glob
import os
import xlwt
from scipy import optimize

def fitFlux(Flux, fitParams):
    A = {}
    B = {}
    fit = {}
    for key in Flux.keys():
        A[key] = []
        for p in fitParams:
            A[key].append(p[key][Flux[key] > 0.2])
        A[key].append(numpy.ones(len(A[key][-1])))
        A[key] = numpy.array(A[key])
        B[key] = numpy.linalg.pinv(A[key])
        fit[key] = B[key].T.dot(Flux[key][Flux[key] > 0.2])
    #print asdf

    return fit
    

fig1 = pyplot.figure(1)
fig1.clear()
ax11 = fig1.add_axes([0.1, 0.1, 0.4, 0.4])
ax11.clear()
ax12 = fig1.add_axes([0.1, 0.5, 0.4, 0.4])
ax12.clear()
ax13 = fig1.add_axes([0.5, 0.1, 0.4, 0.4])
ax13.clear()
ax14 = fig1.add_axes([0.5, 0.5, 0.4, 0.4])
ax14.clear()

fig2 = pyplot.figure(2)
fig2.clear()
ax21 = fig2.add_axes([0.1, 0.1, 0.4, 0.4])
ax21.clear()
ax22 = fig2.add_axes([0.1, 0.5, 0.4, 0.4])
ax22.clear()
ax23 = fig2.add_axes([0.5, 0.1, 0.4, 0.4])
ax23.clear()
ax24 = fig2.add_axes([0.5, 0.5, 0.4, 0.4])
ax24.clear()

fig3 = pyplot.figure(3)
fig3.clear()
ax31 = fig3.add_axes([0.1, 0.1, 0.4, 0.4])
ax31.clear()
ax32 = fig3.add_axes([0.1, 0.5, 0.4, 0.4])
ax32.clear()
ax33 = fig3.add_axes([0.5, 0.1, 0.4, 0.4])
ax33.clear()
ax34 = fig3.add_axes([0.5, 0.5, 0.4, 0.4])
ax34.clear()

fig4 = pyplot.figure(4)
fig4.clear()
ax41 = fig4.add_axes([0.1, 0.1, 0.4, 0.4])
ax41.clear()
ax42 = fig4.add_axes([0.1, 0.5, 0.4, 0.4])
ax42.clear()
ax43 = fig4.add_axes([0.5, 0.1, 0.4, 0.4])
ax43.clear()
ax44 = fig4.add_axes([0.5, 0.5, 0.4, 0.4])
ax44.clear()

#db = CIAO_DatabaseTools.CIAO_Database()

#UTS = [1, 2, 3, 4]
#colors = ['b', 'g', 'r', 'c']

#values = db.query(keywords = ["SEEING", "STREHL", "ASM_SEEING"], AVC_State='ON', timeOfDay='NIGHT', UTS=UTS,
#                  startTime='2017-03-20 00:00:00', endTime='2017-04-01 00:00:00')


FT_Flux = {}
SC_Flux = {}
Strehl = {}
Seeing = {}

FT_Flux[1] = []
FT_Flux[2] = []
FT_Flux[3] = []
FT_Flux[4] = []

SC_Flux[1] = []
SC_Flux[2] = []
SC_Flux[3] = []
SC_Flux[4] = []

Strehl[1] = []
Strehl[2] = []
Strehl[3] = []
Strehl[4] = []

Seeing[1] = []
Seeing[2] = []
Seeing[3] = []
Seeing[4] = []

focus = {}
focus[1] = []
focus[2] = []
focus[3] = []
focus[4] = []

obsTime = []
startTime = {}
startTime[1] = []
startTime[2] = []
startTime[3] = []
startTime[4] = []
datadir = os.environ.get('GRAVITY_DATA')

gravitons = []

for dp, dn, fn in glob.os.walk(datadir):
    for f in fn:
        if 'dualscip2vmred' in f:
            fileBase = dp[len(datadir):]+'/'+f.split('_')[0]
            grav = Graffity.GRAVITY_Data(fileBase=fileBase)
            t = grav.startTime
            obsTime.append(t)

            good = False


ax11.plot(grav.DualSciP2VM.PUPIL_X[1][1::2], grav.DualSciP2VM.PUPIL_Y[1][1::2], color = 'r')
ax11.scatter([0.0], [0.0], color = 'k')
ax11.text(-0.2, 1.2, "GV1")
ax11.set_xbound(-2.0, 2.0)
ax11.set_ybound(-2.0, 2.0)
ax12.plot(grav.DualSciP2VM.PUPIL_X[2][1::2], grav.DualSciP2VM.PUPIL_Y[2][1::2], color = 'r')
ax12.scatter([0.0], [0.0], color = 'k')
ax12.text(-0.2, 1.2, "GV2")
ax12.set_xbound(-2.0, 2.0)
ax12.set_ybound(-2.0, 2.0)
ax13.plot(grav.DualSciP2VM.PUPIL_X[3][1::2], grav.DualSciP2VM.PUPIL_Y[3][1::2], color = 'r')
ax13.scatter([0.0], [0.0], color = 'k')
ax13.text(-0.2, 1.2, "GV3")
ax13.set_xbound(-2.0, 2.0)
ax13.set_ybound(-2.0, 2.0)
ax14.plot(grav.DualSciP2VM.PUPIL_X[4][1::2], grav.DualSciP2VM.PUPIL_Y[4][1::2], color = 'r')
ax14.scatter([0.0], [0.0], color = 'k')
ax14.text(-0.2, 1.2, "GV4")
ax14.set_xbound(-2.0, 2.0)
ax14.set_ybound(-2.0, 2.0)
fig1.show()

ax21.plot(grav.DualSciP2VM.TIME[1::2], grav.DualSciP2VM.PUPIL_X[1][1::2], color = 'b')
ax21.plot(grav.DualSciP2VM.TIME[1::2], grav.DualSciP2VM.PUPIL_Y[1][1::2], color = 'g')
ax21.plot(grav.DualSciP2VM.TIME[1::2], grav.DualSciP2VM.PUPIL_Z[1][1::2], color = 'r')
ax22.plot(grav.DualSciP2VM.TIME[1::2], grav.DualSciP2VM.PUPIL_X[2][1::2], color = 'b')
ax22.plot(grav.DualSciP2VM.TIME[1::2], grav.DualSciP2VM.PUPIL_Y[2][1::2], color = 'g')
ax22.plot(grav.DualSciP2VM.TIME[1::2], grav.DualSciP2VM.PUPIL_Z[2][1::2], color = 'r')
ax23.plot(grav.DualSciP2VM.TIME[1::2], grav.DualSciP2VM.PUPIL_X[3][1::2], color = 'b')
ax23.plot(grav.DualSciP2VM.TIME[1::2], grav.DualSciP2VM.PUPIL_Y[3][1::2], color = 'g')
ax23.plot(grav.DualSciP2VM.TIME[1::2], grav.DualSciP2VM.PUPIL_Z[3][1::2], color = 'r')
ax24.plot(grav.DualSciP2VM.TIME[1::2], grav.DualSciP2VM.PUPIL_X[4][1::2], color = 'b')
ax24.plot(grav.DualSciP2VM.TIME[1::2], grav.DualSciP2VM.PUPIL_Y[4][1::2], color = 'g')
ax24.plot(grav.DualSciP2VM.TIME[1::2], grav.DualSciP2VM.PUPIL_Z[4][1::2], color = 'r')
fig2.show()

ax41.plot(grav.DualSciP2VM.TIME[1::2], grav.DualSciP2VM.PUPIL_U[1][1::2], color = 'b')
ax41.plot(grav.DualSciP2VM.TIME[1::2], grav.DualSciP2VM.PUPIL_V[1][1::2], color = 'g')
ax41.plot(grav.DualSciP2VM.TIME[1::2], grav.DualSciP2VM.PUPIL_W[1][1::2], color = 'r')
ax42.plot(grav.DualSciP2VM.TIME[1::2], grav.DualSciP2VM.PUPIL_U[2][1::2], color = 'b')
ax42.plot(grav.DualSciP2VM.TIME[1::2], grav.DualSciP2VM.PUPIL_V[2][1::2], color = 'g')
ax42.plot(grav.DualSciP2VM.TIME[1::2], grav.DualSciP2VM.PUPIL_W[2][1::2], color = 'r')
ax43.plot(grav.DualSciP2VM.TIME[1::2], grav.DualSciP2VM.PUPIL_U[3][1::2], color = 'b')
ax43.plot(grav.DualSciP2VM.TIME[1::2], grav.DualSciP2VM.PUPIL_V[3][1::2], color = 'g')
ax43.plot(grav.DualSciP2VM.TIME[1::2], grav.DualSciP2VM.PUPIL_W[3][1::2], color = 'r')
ax44.plot(grav.DualSciP2VM.TIME[1::2], grav.DualSciP2VM.PUPIL_U[4][1::2], color = 'b')
ax44.plot(grav.DualSciP2VM.TIME[1::2], grav.DualSciP2VM.PUPIL_V[4][1::2], color = 'g')
ax44.plot(grav.DualSciP2VM.TIME[1::2], grav.DualSciP2VM.PUPIL_W[4][1::2], color = 'r')
fig4.show()

ax31.plot(grav.DITTimes[1][0], grav.FT_Fluxes[1][0],color = 'b')
ax32.plot(grav.DITTimes[2][0], grav.FT_Fluxes[2][0],color = 'b')
ax33.plot(grav.DITTimes[3][0], grav.FT_Fluxes[3][0],color = 'b')
ax34.plot(grav.DITTimes[4][0], grav.FT_Fluxes[4][0],color = 'b')

fig3.show()
