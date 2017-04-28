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
    

fig = pyplot.figure(0)
fig.clear()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
ax.clear()
ax1 = ax.twinx()
ax1.clear()

fig2 = pyplot.figure(1)
fig2.clear()
ax2 = fig2.add_axes([0.1, 0.1, 0.8, 0.8])
ax2.clear()

fig3d = pyplot.figure(2)
fig3d.clear()
ax3d = fig3d.add_axes([0.1, 0.1, 0.8, 0.8], projection='3d')
ax3d.set_xlabel("Seeing")
ax3d.set_ylabel("Strehl")
ax3d.set_zlabel("Flux")

db = CIAO_DatabaseTools.CIAO_Database()

UTS = [1, 2, 3, 4]
colors = ['b', 'g', 'r', 'c']

values = db.query(keywords = ["SEEING", "STREHL", "ASM_SEEING"], AVC_State='ON', timeOfDay='NIGHT', UTS=UTS,
                  startTime='2017-03-20 00:00:00', endTime='2017-04-01 00:00:00')


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
        if 'astroreduced' in f:
            fileBase = dp[len(datadir):]+'/'+f.split('_')[0]
            grav = Graffity.GRAVITY_Data(fileBase=fileBase)
            t = grav.startTime
            obsTime.append(t)

            good = False

            order = numpy.argsort(numpy.abs(values[1][:,-4] - t))
            if numpy.abs(t- values[1][order[0],-4]) < 300.0:
                good = True
                grav.addAOData(UT=1, Strehl=values[1][order[0],1], Seeing=values[1][order[0],0], 
                               ASM_Seeing=values[1][order[0],2])
                FT_Flux[1].append(numpy.mean(grav.FT_Fluxes[1]))
                SC_Flux[1].append(numpy.mean(grav.SC_Fluxes[1]))
                Strehl[1].append(values[1][order[0],1])
                Seeing[1].append(values[1][order[0],0])
                DL = Graffity.DataLogger(directory=values[1][order[0],-3])
                DL.loadData()
                focus[1].append(DL.getRefSlopeZernikes()[2][0])
                startTime[1].append(grav.startTime)

            order = numpy.argsort(numpy.abs(values[2][:,-4] - t))
            if numpy.abs(t- values[2][order[0],-4]) < 300.0:
                good = True
                grav.addAOData(UT=2, Strehl=values[2][order[0],1], Seeing=values[2][order[0], 0],
                               ASM_Seeing=values[2][order[0],2])
                FT_Flux[2].append(numpy.mean(grav.FT_Fluxes[2]))
                SC_Flux[2].append(numpy.mean(grav.SC_Fluxes[2]))
                Strehl[2].append(values[2][order[0],1])
                Seeing[2].append(values[2][order[0],0])
                DL = Graffity.DataLogger(directory=values[2][order[0],-3])
                DL.loadData()
                focus[2].append(DL.getRefSlopeZernikes()[2][0])
                startTime[2].append(grav.startTime)


            order = numpy.argsort(numpy.abs(values[3][:,-4] - t))
            if numpy.abs(t- values[3][order[0],-4]) < 300.0:
                good = True
                grav.addAOData(UT=3, Strehl=values[3][order[0],1], Seeing=values[3][order[0], 0],
                               ASM_Seeing=values[3][order[0],2])
                FT_Flux[3].append(numpy.mean(grav.FT_Fluxes[3]))
                SC_Flux[3].append(numpy.mean(grav.SC_Fluxes[3]))
                Strehl[3].append(values[3][order[0],1])
                Seeing[3].append(values[3][order[0],0])
                DL = Graffity.DataLogger(directory=values[3][order[0],-3])
                DL.loadData()
                focus[3].append(DL.getRefSlopeZernikes()[2][0])
                startTime[3].append(grav.startTime)


            order = numpy.argsort(numpy.abs(values[4][:,-4] - t))
            if numpy.abs(t- values[4][order[0],-4]) < 300.0:
                good = True
                grav.addAOData(UT=4, Strehl=values[4][order[0],1], Seeing=values[4][order[0], 0],
                               ASM_Seeing=values[4][order[0],2])
                FT_Flux[4].append(numpy.mean(grav.FT_Fluxes[4]))
                SC_Flux[4].append(numpy.mean(grav.SC_Fluxes[4]))
                Strehl[4].append(values[4][order[0],1])
                Seeing[4].append(values[4][order[0],0])
                DL = Graffity.DataLogger(directory=values[4][order[0],-3])
                DL.loadData()
                focus[4].append(DL.getRefSlopeZernikes()[2][0])
                startTime[4].append(grav.startTime)

            if good:
                gravitons.append(grav)


obsTime = numpy.array(obsTime)
obsTime -= numpy.min(obsTime)
obsTime /= numpy.max(obsTime)
obsTime *= 35.0
obsTime += 15.0

FT_Flux[1] = numpy.array(FT_Flux[1])
FT_Flux[1] /= numpy.max(FT_Flux[1])
FT_Flux[2] = numpy.array(FT_Flux[2])
FT_Flux[2] /= numpy.max(FT_Flux[2])
FT_Flux[3] = numpy.array(FT_Flux[3])
FT_Flux[3] /= numpy.max(FT_Flux[3])
FT_Flux[4] = numpy.array(FT_Flux[4])
FT_Flux[4] /= numpy.max(FT_Flux[4])
SC_Flux[1] = numpy.array(SC_Flux[1])
SC_Flux[1] /= numpy.max(SC_Flux[1])
SC_Flux[2] = numpy.array(SC_Flux[2])
SC_Flux[2] /= numpy.max(SC_Flux[2])
SC_Flux[3] = numpy.array(SC_Flux[3])
SC_Flux[3] /= numpy.max(SC_Flux[3])
SC_Flux[4] = numpy.array(SC_Flux[4])
SC_Flux[4] /= numpy.max(SC_Flux[4])
Strehl[1] = numpy.array(Strehl[1])
Strehl[2] = numpy.array(Strehl[2])
Strehl[3] = numpy.array(Strehl[3])
Strehl[4] = numpy.array(Strehl[4])
Seeing[1] = numpy.array(Seeing[1])
Seeing[2] = numpy.array(Seeing[2])
Seeing[3] = numpy.array(Seeing[3])
Seeing[4] = numpy.array(Seeing[4])

#fit = fitFlux(FT_Flux, [Strehl, Seeing])
#fit = fitFlux(SC_Flux, [Strehl])
fit = fitFlux(FT_Flux, [Strehl])
print fit

colors = ['b', 'g', 'r', 'y']
for UT in UTS:
    startTime[UT] = numpy.array(startTime[UT])
    startTime[UT] -= numpy.min(startTime[UT])
    startTime[UT] /= 3600.0
    FT_Flux[UT] = numpy.array(FT_Flux[UT])
    SC_Flux[UT] = numpy.array(SC_Flux[UT])
    Strehl[UT] = numpy.array(Strehl[UT])
    Seeing[UT] = numpy.array(Seeing[UT])
    #FitFlux = fit[UT].dot([Strehl[UT], Seeing[UT], numpy.ones(len(Seeing[UT]))])
    FitFlux = fit[UT].dot([Strehl[UT], numpy.ones(len(Strehl[UT]))])
    order = numpy.argsort(startTime[UT])
    ax.scatter(startTime[UT], FT_Flux[UT], c=colors[UT-1])
    #ax.scatter(startTime[UT], SC_Flux[UT], c=colors[UT-1])
    ax.plot(startTime[UT][order], FT_Flux[UT][order], c=colors[UT-1], lw=2.0, ls='-')
    #ax.plot(startTime[UT][order], SC_Flux[UT][order], c=colors[UT-1], lw=2.0, ls='-')
    ax1.plot(startTime[UT][order], Strehl[UT][order], c=colors[UT-1], lw=1.0, ls='--')
    ax1.plot(startTime[UT][order], Seeing[UT][order], c=colors[UT-1], lw=1.0, ls='-.')
    ax.plot(startTime[UT][order], FitFlux[order], c=colors[UT-1], lw=3.0, ls=':')
    ax3d.scatter(xs=Seeing[UT], ys=Strehl[UT], zs=FT_Flux[UT], c=colors[UT-1])
    #ax3d.scatter(xs=Seeing[UT], ys=Strehl[UT], zs=SC_Flux[UT], c=colors[UT-1])
    ax2.scatter(Strehl[UT], FT_Flux[UT], c=colors[UT-1])
    #ax2.scatter(Strehl[UT], SC_Flux[UT], c=colors[UT-1])
    ax2.scatter(Strehl[UT], FitFlux, c=colors[UT-1], marker = "x")

ax.set_xlabel("Strehl Ratio")
ax.set_ylabel("FT Flux")
#ax.set_xbound(1.0, 0.0)
#ax.set_ybound(0, 1.8e6)

#raw_input()

#grav.plotReducedFTFluxes(ax=ax)

fig.show()
fig3d.show()
fig2.show()



"""
Book = xlwt.Workbook()
sheet = Book.add_sheet("Strehl")
sheet.write(0, 0, "Archive File")
sheet.write(0, 1, "Strehl UT1")
sheet.write(0, 2, "Strehl UT2")
sheet.write(0, 3, "Strehl UT3")
sheet.write(0, 4, "Strehl UT4")
i = 1
for grav in gravitons:
    sheet.write(i, 0, grav.AstroReduced.masterheader.get('ARCFILE'))
    if 1 in grav.Strehl.keys():
        sheet.write(i, 1, grav.Strehl[1])
    else:
        sheet.write(i, 1, numpy.nan)
    if 2 in grav.Strehl.keys():
        sheet.write(i, 2, grav.Strehl[2])
    else:
        sheet.write(i, 2, numpy.nan)
    if 3 in grav.Strehl.keys():
        sheet.write(i, 3, grav.Strehl[3])
    else:
        sheet.write(i, 3, numpy.nan)
    if 4 in grav.Strehl.keys():
        sheet.write(i, 4, grav.Strehl[4])
    else:
        sheet.write(i, 4, numpy.nan)
    i += 1

Book.save('Strehl.xls')
"""
