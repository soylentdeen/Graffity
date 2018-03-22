import astropy.io.fits as pyfits
import numpy
import sys
sys.path.append('../')

import Graffity
import PlotTools
import os
import os.path

strehls = {}


Z2S = pyfits.getdata('cimdatZernike2Slopes.fits')
S2Z = numpy.linalg.pinv(Z2S)

PlotTools.clearAllPlots()
figs, axes = PlotTools.configurePlots(3)

residuals = {}
tiptilt = {}
for i in [1, 2, 3, 4]:
    residuals[i] = []
    tiptilt[i] = []
    strehls[i] = []
    for base, dirs, files in os.walk('/tera/5/CIAO/CIAO%d/' % i):
        if ('CIAO_LOOP_0001.fits' in files) and ('DATA_EXPO' in base):
            print base
            head = pyfits.getheader(os.path.join(base, 'CIAO_LOOP_0001.fits'))
            gradients = pyfits.getdata(os.path.join(base,'CIAO_LOOP_0001.fits')).field('Gradients')
            strehls[i].append(head.get('ESO AOS ATM SR'))
            resids = []
            tt = []
            for g in gradients:
                resids.append(S2Z.T.dot(g))
                tt.append([numpy.mean(g[0::2]), numpy.mean(g[1::2])])
            residuals[i].append(numpy.mean(numpy.array(resids)**2.0, axis=0)**0.5)
            tiptilt[i].append(numpy.std(tt, axis=0)*0.5)

    tiptilt[i] = numpy.array(tiptilt[i])
    residuals[i] = numpy.array(residuals[i])
    strehls[i] = numpy.array(strehls[i])

color = {1:'b', 2:'g', 3:'r', 4:'y'}

for i in [1, 2, 3, 4]:
    axes[0].scatter(strehls[i], residuals[i][:,0], color = color[i], label =
            "UT%d" % i)
    axes[0].scatter(strehls[i], residuals[i][:,1], color = color[i])
#axes[1].scatter(strehls, residuals[:,2], color = 'b')
#axes[1].scatter(strehls, residuals[:,3], color = 'g')
#axes[1].scatter(strehls, residuals[:,4], color = 'r')
#axes[1].scatter(strehls, residuals[:,5], color = 'y')
#axes[2].scatter(strehls, tiptilt[:,0], color = 'b')
#axes[2].scatter(strehls, tiptilt[:,1], color = 'r')
axes[0].set_xlabel("Kband Strehl Ratio")
axes[0].set_ylabel("Tip/Tilt Residuals (arcsec)")
axes[0].set_ybound(0.0, 0.02)
axes[0].set_xbound(0.0, 1.0)
axes[0].legend(loc=1, scatterpoints=1)
figs[0].suptitle("Tip Tilt Residuals from CIAO Data")
figs[0].show()
figs[0].savefig("TipTilt_Residuals.png")
#figs[1].show()
#figs[2].show()
