import sys
sys.path.append('../')
import Graffity
import numpy
import PlotTools

figs, axes = PlotTools.configurePlots(1)

filebase = '/media/cdeen/My Passport/GRAVITY/GRAVITY/reduced_20180307/GRAVI.2017-07-10T00:40:46.319'

Grav = Graffity.GRAVITY_Dual_Sci_P2VM(fileBase=filebase, processAcqCamData=True)


color = {0:'b', 1:'g', 2:'r', 3:'y'}
for i in Grav.AcqCamDat.newStrehl.data.keys():
    #axes[0].scatter(Grav.AcqCamDat.newStrehl.data[i], Grav.AcqCamDat.TOTALFLUX_SC.data[i], color=color[i])
    angle = numpy.deg2rad(Grav.NORTH_ANGLE[i])
    xy = numpy.array([Grav.AcqCamDat.newSC_FIBER_DX.data[i], Grav.AcqCamDat.newSC_FIBER_DY.data[i]])*Grav.AcqCamDat.newSCALE.data[i]
    rotmat = numpy.array([[numpy.cos(angle), -numpy.sin(angle)],[numpy.sin(angle), numpy.cos(angle)]])
    radec = rotmat.dot(xy)
    radec -= numpy.array([numpy.mean(radec, axis=1)]).T
    radec += numpy.array([[Grav.MET_SOBJ_DRA[i]],[Grav.MET_SOBJ_DDEC[i]]])
    FiberOffset = numpy.sum(radec**2.0, axis=0)**0.5
    axes[0].scatter(FiberOffset, Grav.AcqCamDat.TOTALFLUX_SC.data[i], color=color[i])
    #axes[0].scatter(numpy.mean(Grav.AcqCamDat.SC_FIBER_DX.data[i]*Grav.AcqCamDat.SCALE.data[i]), numpy.mean(Grav.AcqCamDat.SC_FIBER_DY.data[i]*Grav.AcqCamDat.SCALE.data[i]), color=color[i])

figs[0].show()
