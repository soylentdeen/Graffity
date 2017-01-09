import scipy
import numpy
import scipy.io as sio
import astropy.io.fits as pyfits

import matplotlib.pyplot as pyplot

fig = pyplot.figure(0)
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

datafile ='/home/cdeen/Data/CIAO/UT3/TemplateData/2016-08-31_3/ANAMORPH-205133/IMs/Demodulated/IMs.mat'
DemodulatedIMs = sio.loadmat(datafile, struct_as_record=False, squeeze_me=True)

TT2HO = pyfits.getdata('/home/cdeen/Code/Matlab/Database/UT3/CDMS/Common/RecnOptimiser.TT2HO.fits')

iTT2HO = numpy.linalg.linalg.pinv(TT2HO)

IMs = DemodulatedIMs['IMs']
Headers = DemodulatedIMs['DataHeader']

nAct = 60
nSubAp = 68
nAx = 2
nAngles = IMs.shape[-1]

HOIMs = IMs[:,0:60,:].swapaxes(1,2)


HOIMs = HOIMs.reshape([nSubAp*nAx*nAngles, nAct])



TTRemover = numpy.lib.twodim_base.eye(nAct) - TT2HO.dot(iTT2HO)

TTFHOIMs = HOIMs.dot(TTRemover)

u, s, Basis = numpy.linalg.svd(TTFHOIMs)

Basis = Basis.T

HOBasis = Basis[:58]
TTBasis = Basis[58:]

HOCoefficients = iTT2HO.dot(HOBasis.T)
TTCoefficients = iTT2HO.dot(TTBasis.T)
Contamination = norm(HOCoefficients)/norm(TTCoefficients)
if Contamination > 1e-5:
    print('Failed to extract the TT Modes');


ax.matshow(HOBasis)

fig.show()
