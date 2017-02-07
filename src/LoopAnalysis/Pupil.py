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

fig0 = pyplot.figure(0)
fig0.clear()
ax0 = fig0.add_axes([0.1, 0.1, 0.8, 0.8])

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

for ciao in [1, 2, 3, 4]:
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
            vcmu.append(CB.header.get('ESO STS VCM2 POSX'))
            vcmw.append(CB.header.get('ESO STS VCM2 POSY'))
    Time[ciao] = numpy.array(t)
    Flux[ciao] = numpy.array(flux) - numpy.array([1.0, 1.0])
    Alt[ciao] = numpy.array(alt)
    Az[ciao] = numpy.array(az)
    VCMU[ciao] = numpy.array(vcmu)
    VCMW[ciao] = numpy.array(vcmw)


for ciao in [1, 2, 3, 4]:
    fit, matrix = doFit(Flux[ciao], Alt[ciao], Az[ciao], VCMU[ciao], VCMW[ciao])
    ax0.clear()
    ax0.scatter(Flux[ciao][:,0], Flux[ciao][:,1], color = 'b')
    ax0.scatter(fit[0,:], fit[1,:], color = 'g')
    fig0.show()
    input()

#for ciao in [1, 2, 3, 4]:
#    ax0.clear()
#    ax0.scatter(Flux[ciao][:,0], Flux[ciao][:,1], s=Az[ciao], c=Alt[ciao])
#    fig0.show()
#    fig0.savefig("PupilCentroid_"+str(ciao)+".png")
#    input()
