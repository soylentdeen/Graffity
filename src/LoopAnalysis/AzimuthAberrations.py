import Graffity
import numpy
import matplotlib.pyplot as pyplot
import scipy
import pyfits
import sys
import glob
import time
from scipy import optimize

datadir = sys.argv[1]

pyplot.rcParams['font.size'] = 16.0
pyplot.rcParams['legend.fontsize'] = 10.0
fig = pyplot.figure(0, figsize = (12, 8), dpi=300)
fig.clear()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

dataloggers = glob.glob(datadir+'DATA_LOGGER*')

def fitSinusoid(derot, amplitude):
    #errfunc = lambda p, th: p[0] + p[1]*numpy.cos(p[2]*th*3.14159/180.0 + p[3])
    errfunc = lambda p, th: p[0] + p[1]*numpy.cos(4*th*3.14159/180.0 + p[2])
    fitfunc = lambda p, y, th: numpy.abs(errfunc(p, th) - y)
    guess = [numpy.average(amplitude), numpy.std(amplitude), 0.0]
    pfit = optimize.leastsq(fitfunc, guess, args = (amplitude, derot))
    angles = numpy.arange(180)
    fit = errfunc(pfit[0], angles)
    return angles, fit, pfit[0][1], pfit[0][2], pfit[0][2]

def getSaturationRate(frames):
    satFrames = []
    for frame in frames:
        satFrames.append(float(len(frame[abs(frame) > 0.35]))/float(len(frame)))

    return numpy.mean(satFrames)

TTM = []
Azimuth = []
Derot = []
VCM_U = []
VCM_W = []
TTM_RMS = []
HODM_SHAPE = []
HODM_ZERN = []
Slope_RMS = []
Saturations = []

Noll = ['Piston', 'Tip', 'Tilt', 'Defocus', 'Oblique Ast.', 'Vert. Ast.', 'Vert. Coma',
        'Horiz. Coma', 'Vert. Tref.', 'Oblique Tref.', 'Spherical']
Selected = {}
Selected['Piston']=False
Selected['Tip'] = True
Selected['Tilt'] = True
Selected['Defocus'] = False
Selected['Oblique Ast.'] = True
Selected['Vert. Ast.'] = True
Selected['Vert. Coma'] = False
Selected['Horiz. Coma'] = False
Selected['Vert. Tref.'] = False
Selected['Oblique Tref.'] = False
Selected['Spherical'] = False

dataloggers.sort()

Az = -90

templateTime = []

for datalogger in dataloggers[2:-3]:
    if len(Azimuth) == 0:
        Azimuth.append(Az)
    elif len(Azimuth) == 10:
        Azimuth.append(Azimuth[-1]+20)
    else:
        Azimuth.append(Azimuth[-1]+10)
    CB = Graffity.CircularBuffer(df=datalogger+'/CIAO_LOOP_0001.fits', 
            Z2DM='../../data/cimdatZernike2DM.fits')

    CB.loadTemplateMaps()

    Derot.append(CB.header.get("HIERARCH ESO INS DROT POSANG"))
    templateTime.append(time.mktime(time.strptime(CB.header.get(
        "HIERARCH ESO TPL START"), '%Y-%m-%dT%H:%M:%S')))
    TTM.append(numpy.average(CB.TTM, axis=0))

    VCM_U.append(CB.header.get("HIERARCH ESO STS VCM2 GUIDE U"))
    VCM_W.append(CB.header.get("HIERARCH ESO STS VCM2 GUIDE W"))

    HODM_SHAPE.append(numpy.average(CB.HODM, axis=0) - CB.HO_ACT_REF_MAP[0])
    HODM_ZERN.append(CB.calculateMirrorZernikes(HODM_SHAPE[-1])/1e-6)

    TTM_RMS.append(numpy.std(CB.TTM, axis=0))

    Slope_RMS.append(numpy.std(CB.Gradients, axis=0))

    Saturations.append(getSaturationRate(CB.HODM))

mirrorFlat = CB.calculateMirrorZernikes(CB.HO_ACT_REF_MAP[0])/1e-6

templateTime = numpy.array(templateTime)
templateTime = (templateTime - numpy.min(templateTime))/60.0
TTM = numpy.array(TTM)
TTM_RMS = numpy.array(TTM_RMS)
HODM_SHAPE = numpy.array(HODM_SHAPE)
HODM_ZERN = numpy.array(HODM_ZERN)
Slope_RMS = numpy.array(Slope_RMS)
Saturations = numpy.array(Saturations)

# Tip Tilt Position

# Field Lens Positions

# VCM Position

# TTM RMS

# Slope RMS

# Saturations


#for HODM in HODM_ZERN:
#    ax.plot(HODM)

Derot = numpy.array(Derot)
order = numpy.argsort(Derot)

breakTime = 20.0

early = (templateTime - breakTime) < 0.0
late = (templateTime - breakTime) > 0.0

for i in range(len(Noll)):
    if Selected[Noll[i]]:
        angles, fit, amplitude, omega, phi = fitSinusoid(Derot[late], HODM_ZERN[late,i])
        blah = ax.plot(Derot[late], HODM_ZERN[late,i], label = Noll[i], lw=3.0)
                #lw=2.0, color = numpy.random.rand(3,1))
        ax.plot(angles, fit, color = blah[0].get_c(), ls = '--', lw=3.0)
        print("%s %.2f %.2f\n" % (Noll[i], amplitude, omega*180.0/3.14159))

ax.legend(ncol=6)
#ax.scatter(TTM[:,0], TTM[:,1])
#ax.scatter(Azimuth, TTM[:,0])
#ax.scatter(Azimuth, TTM[:,1])
#ax.scatter(Azimuth, numpy.average(TTM_RMS, axis=1))
#ax.scatter(Azimuth, numpy.std(Slope_RMS, axis=1))
#ax.scatter(Azimuth, Saturations)
#ax.set_xbound(-0.1, 0.1)
#ax.set_ybound(-0.1, 0.1)
#ax.set_aspect('equal')
ax.set_title("Mirror Aberrations under Derotator Rotation")
ax.set_xlabel("Derotator position (Degrees)")
ax.set_ylabel(r"Zernike Coefficient ($\mu$m RMS)")
fig.savefig("Derotator_v_MirrorZern.png")

ax.clear()
ax.bar(numpy.arange(len(mirrorFlat))+1, mirrorFlat)
ax.set_title("CIAO Flat Pattern Shape")
ax.set_xlabel("Noll Index")
ax.set_ylabel(r"Zernike Coefficient ($\mu$m RMS)")
fig.show()
fig.savefig("CIAO_FLAT_DM.png")
