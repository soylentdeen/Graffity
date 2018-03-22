import sys
sys.path.append('../')
import Graffity
import numpy
import PlotTools
import glob
import astropy.io.fits as pyfits
from scipy import optimize
from os import walk
from os.path import join

def fitGaussian(x, y):
    errfunc = lambda p, x, y: numpy.abs(p[0]*numpy.exp(-(x**2.0)/(2*p[1]**2.0)) - y)
    coeffs = [300000.0, 25.0]
    pfit, success = optimize.leastsq(errfunc, coeffs, args=(x,y))
    return pfit


figs, axes = PlotTools.configurePlots(4)

observationFiles = numpy.array([])
startdir = '/gvstore1/forFrank/'

for root, dirs, files in walk(startdir):
    if 'reduced_20180307' in root:
        for f in files:
            if 'dualscip2vmred.fits' in f:
                observationFiles = numpy.append(observationFiles, join(root, f))

color = {0:'b', 1:'g', 2:'r', 3:'y'}

offsets = {}
offsets[0] = numpy.array([])
offsets[1] = numpy.array([])
offsets[2] = numpy.array([])
offsets[3] = numpy.array([])
strehls = {}
strehls[0] = numpy.array([])
strehls[1] = numpy.array([])
strehls[2] = numpy.array([])
strehls[3] = numpy.array([])
flux = {}
flux[0] = numpy.array([])
flux[1] = numpy.array([])
flux[2] = numpy.array([])
flux[3] = numpy.array([])
for f in observationFiles:
    SObjName = pyfits.getheader(f).get("ESO INS SOBJ NAME")
    if SObjName == 'IRS16NW':
        Grav = Graffity.GRAVITY_Dual_Sci_P2VM(fileBase=f[:-20], processAcqCamData=True)
        print "Good - %s" %f
        for i in Grav.AcqCamDat.newStrehl.data.keys():
            FiberOffset = (Grav.MET_SOBJ_DRA[i]**2.0 + Grav.MET_SOBJ_DDEC[i]**2.0)**0.5
            offsets[i] = numpy.append(offsets[i], FiberOffset)
            strehls[i] = numpy.append(strehls[i], Grav.Strehl[i])
            flux[i] = numpy.append(flux[i], numpy.mean(Grav.AcqCamDat.TOTALFLUX_SC.data[i]))
            axes[2].scatter([FiberOffset],[flux[i][-1]],
                    color=color[i], s=2*2**(10*strehls[i][-1]))
            axes[3].scatter([strehls[i][-1]], [flux[i][-1]], color=color[i])

        del(Grav)
    else:
        print "Bad"

binSize=0.2
StrehlBins = numpy.arange(0.05, 0.65, binSize)
offsetRange = numpy.linspace(0.0, 50.0)

Amps = {}
sigs = {}
srs = {}

dFib = numpy.linspace(0, 50.0)

for i in color.keys():
    Amps[i] = []
    sigs[i] = []
    srs[i] = []
    for s in StrehlBins:
        current = (strehls[i] > s) & (strehls[i] < s+binSize)
        if numpy.sum(current) > 5:
            #axes[0].clear()
            #axes[0].scatter(offsets[i][current], flux[i][current], color=color[i])
            fit = fitGaussian(offsets[i][current], flux[i][current])
            #print fit[1]
            #axes[0].plot(offsetRange, fit[0]*numpy.exp(-(offsetRange)**2.0/(2.0*fit[1])**2.0))
            #figs[0].show()
            #raw_input()
            Amps[i].append(fit[0])
            sigs[i].append(fit[1])
            srs[i].append(s)
            if sigs[i][-1] < 70:
                axes[2].plot(dFib,
                        Amps[i][-1]*numpy.exp(-(dFib)**2.0/(2.0*sigs[i][-1])**2.0),
                        color = color[i], lw=s)

    Amps[i] = numpy.array(Amps[i])
    sigs[i] = numpy.array(sigs[i])
    srs[i] = numpy.array(srs[i])

    axes[0].scatter(srs[i], Amps[i]/1000000.0, color=color[i], label="Telescope = %d"%(i+1))
    axes[1].scatter(srs[i][sigs[i] < 70], sigs[i][sigs[i] < 70],
            color=color[i], label="Telescope = %d"%(i+1))

axes[0].set_xlabel("Strehl Ratio")
axes[0].set_ylabel("Amplitude")
axes[1].set_xlabel("Strehl Ratio")
axes[1].set_ylabel("Sigma (mas)")
axes[2].set_xlabel("Fiber Offset (mas)")
axes[2].set_ylabel("Total Flux (SC)")
axes[3].set_xlabel("Strehl Ratio")
axes[3].set_ylabel("Total Flux (SC)")
figs[0].suptitle("IRS16NW")
figs[1].suptitle("IRS16NW")
figs[2].suptitle("IRS16NW")
figs[3].suptitle("IRS16NW")
axes[0].legend(loc=4)
axes[1].legend(loc=1)
figs[0].show()
figs[0].savefig("Amplitude_V_Strehl.png")
figs[1].show()
figs[1].savefig("Sigma_V_Strehl.png")
figs[2].show()
figs[3].show()
figs[2].savefig("Offset_V_Flux.png")
figs[3].savefig("Strehl_V_Flux.png")

