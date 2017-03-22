import scipy
import numpy
from matplotlib import pyplot
import Graffity
import glob
import time
from scipy import ndimage

def doFit(centroids, alt, az, u, w):
    #Ax = b
    A = numpy.matrix(numpy.array([numpy.ones(len(alt)), alt, az, u, w]).T)
    x = A.I.dot(centroids)
    fitPoints = x.T.dot(A.T)
    return fitPoints, x

def doAltAzFit(centroids, alt, az):
    #Ax = b
    A = numpy.matrix(numpy.array([numpy.ones(len(az)), alt, az]).T)
    x = A.I.dot(centroids)
    fitPoints = x.T.dot(A.T)
    return fitPoints, x

fig0 = pyplot.figure(0)
fig0.clear()
ax01 = fig0.add_axes([0.1, 0.1, 0.4, 0.4])
ax02 = fig0.add_axes([0.1, 0.5, 0.4, 0.4])
ax03 = fig0.add_axes([0.5, 0.1, 0.4, 0.4])
ax04 = fig0.add_axes([0.5, 0.5, 0.4, 0.4])
fig1 = pyplot.figure(1)
fig1.clear()
ax11 = fig1.add_axes([0.1, 0.1, 0.4, 0.4])
ax12 = fig1.add_axes([0.1, 0.5, 0.4, 0.4])
ax13 = fig1.add_axes([0.5, 0.1, 0.4, 0.4])
ax14 = fig1.add_axes([0.5, 0.5, 0.4, 0.4])

CIAO_dir = '/home/cdeen/Data/CIAO/JanComm/2017-01-08_'
CDMS_Dir = '/home/cdeen/Data/CIAO/DATABASE/UT4/OFFAXIS/Matrices/Common/'
CDMS_BaseDir = '/home/cdeen/Data/CIAO/DATABASE/UT'
CDMS_ConfigDir = '/OFFAXIS/Matrices/Common/'

Az = {}
Alt = {}
Flux = {}
Time = {}
VCMU = {}
VCMW = {}

startTime = time.mktime(time.strptime('2016-09-21T23:00:00', '%Y-%m-%dT%H:%M:%S'))

for ciao, ax in zip([1, 2, 3, 4], [ax11, ax12, ax13, ax14]):
    CIAO_datadir = CIAO_dir+str(ciao)+'/'
    az = []
    alt = []
    flux = []
    t = []
    vcmu = []
    vcmw = []
    for dp, dn, fn in glob.os.walk(CIAO_datadir):
        if ((len(dn) == 0) and ('DATA_LOGGER' in dp)):    #26, 34, 35, 43
            CB =Graffity.CircularBuffer(dp+'/CIAO_LOOP_0001.fits', CDMS_BaseDir=CDMS_BaseDir, 
                 CDMS_ConfigDir = CDMS_ConfigDir, S2M=dp+'/RecnOptimiser.S2M_0001.fits', 
                 ModalBasis='RecnOptimiser.ModalBasis.fits',Z2DM='RecnOptimiser.Z2DM.fits', 
                 CM=dp+'/Recn.REC1.CM_0001.fits', HOIM=dp+'/RecnOptimiser.HO_IM_0001.fits',
                 TT2HO='RecnOptimiser.TT2HO.fits', DM2Z='RecnOptimiser.DM2Z.fits',
                 TTM2Z='RecnOptimiser.TTM2Z.fits', prefix=str(ciao))

            tm = (time.mktime(time.strptime(CB.header.get('ESO TPL START'), '%Y-%m-%dT%H:%M:%S'))-startTime)/3600.0
            print("file: %s " % dp)
            image = Graffity.WFS_Frame(intensities = numpy.mean(CB.Intensities, axis=0))
            flux.append(ndimage.measurements.center_of_mass(image.image[3:6,3:6]))
            t.append(tm)
            alt.append(CB.header.get('ESO TEL ALT'))
            az.append(CB.header.get('ESO TEL AZ'))
            vcmu.append(CB.header.get('ESO STS VCM2 GUIDE U'))
            vcmw.append(CB.header.get('ESO STS VCM2 GUIDE W'))
    Time[ciao] = numpy.array(t)
    Flux[ciao] = numpy.array(flux) - numpy.array([1.0, 1.0])
    Alt[ciao] = numpy.array(alt)
    Az[ciao] = numpy.array(az)
    VCMU[ciao] = numpy.array(vcmu)
    VCMW[ciao] = numpy.array(vcmw)
    image.plotIntensities(ax=ax)


for ciao, ax in zip([1, 2, 3, 4], [ax01, ax02, ax03, ax04]):
    fit, matrix = doFit(Flux[ciao], Alt[ciao], Az[ciao], VCMU[ciao], VCMW[ciao])
    #fit, matrix = doAltAzFit(Flux[ciao], Alt[ciao], Az[ciao])
    ax.clear()
    ax.scatter(Flux[ciao][:,0], Flux[ciao][:,1], color = 'b')
    ax.scatter(fit[0,:], fit[1,:], color = 'g')
    ax.text(0.1, 0.9, 'UT%d'%ciao, transform = ax.transAxes)

fig0.suptitle("Pupil Decenter")
fig0.show()

fig1.suptitle("Pupil Illumination")
fig1.savefig("Pupil.png")
fig0.savefig("decenter.png")

#for ciao in [1, 2, 3, 4]:
#    ax0.clear()
#    ax0.scatter(Flux[ciao][:,0], Flux[ciao][:,1], s=Az[ciao], c=Alt[ciao])
#    fig0.show()
#    fig0.savefig("PupilCentroid_"+str(ciao)+".png")
#    input()
