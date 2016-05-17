import Graffity
import numpy
import matplotlib.pyplot as pyplot
import scipy
import pyfits
import sys
import glob

datadir = sys.argv[1]

pyplot.rcParams['font.size'] = 20.0
pyplot.rcParams['legend.fontsize'] = 12.0
fig = pyplot.figure(0, figsize = (8, 6), dpi=300)
fig.clear()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

dataloggers = glob.glob(datadir+'DATA_LOGGER*')

def getSaturationRate(frames):
    satFrames = []
    for frame in frames:
        satFrames.append(float(len(frame[abs(frame) > 0.35]))/float(len(frame)))

    return numpy.mean(satFrames)

TTM = []
Azimuth = []
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
Selected['Defocus'] = True
Selected['Oblique Ast.'] = True
Selected['Vert. Ast.'] = True
Selected['Vert. Coma'] = False
Selected['Horiz. Coma'] = False
Selected['Vert. Tref.'] = False
Selected['Oblique Tref.'] = True
Selected['Spherical'] = False

dataloggers.sort()

Az = -90

for datalogger in dataloggers[2:-3]:
    if len(Azimuth) == 0:
        Azimuth.append(Az)
    elif len(Azimuth) == 10:
        Azimuth.append(Azimuth[-1]+20)
    else:
        Azimuth.append(Azimuth[-1]+10)
    CB = Graffity.CircularBuffer(df=datalogger+'/CIAO_LOOP_0001.fits', 
            Z2DM='../../data/cimdatZernike2DM.fits')

    TTM.append(numpy.average(CB.TTM, axis=0))

    VCM_U.append(CB.header.get("HIERARCH ESO STS VCM2 GUIDE U"))
    VCM_W.append(CB.header.get("HIERARCH ESO STS VCM2 GUIDE W"))

    HODM_SHAPE.append(numpy.average(CB.HODM, axis=0))
    HODM_ZERN.append(CB.calculateMirrorZernikes(HODM_SHAPE[-1])/1e-6)

    TTM_RMS.append(numpy.std(CB.TTM, axis=0))

    Slope_RMS.append(numpy.std(CB.Gradients, axis=0))

    Saturations.append(getSaturationRate(CB.HODM))

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


for i in range(len(Noll)):
    if Selected[Noll[i]]:
        ax.plot(Azimuth, HODM_ZERN[:,i], label = Noll[i], lw=3.0)
                #lw=2.0, color = numpy.random.rand(3,1))

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
ax.set_title("Mirror Aberrations under Azimuth Rotation")
ax.set_xlabel("Azimuth position (Degrees)")
ax.set_ylabel(r"Zernike Coefficient ($\mu$m RMS)")
fig.show()
fig.savefig("Azimuth_v_MirrorZern.png")
