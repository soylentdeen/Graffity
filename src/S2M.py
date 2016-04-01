import scipy
import pyfits
import scipy.linalg
import numpy

z2s = pyfits.getdata("../data/Z2S_136s_119z.fits")

s2z = scipy.linalg.pinv(z2s)

blah = numpy.array(s2z[:,:50], dtype=numpy.float32)

s2Zernike = pyfits.PrimaryHDU(data=blah.T)

s2Zernike.writeto("slopes2zerns.fits", clobber=True)

S2Z = pyfits.getdata("Slopes2Z.fits")

junk = numpy.array(S2Z[:50,:]*1e5, dtype=numpy.float32)

s2Zernike = pyfits.PrimaryHDU(data=junk)

s2Zernike.writeto("S2Z.fits", clobber=True)
