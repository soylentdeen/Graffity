import matplotlib.pyplot as pyplot
import astropy.io.fits as pyfits

fig = pyplot.figure(0)
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])


good = pyfits.getdata('/home/cdeen/Data/CIAO/UT4/GoodDatabase/UT4/CDMS/Common/RecnOptimiser.ModalBasis.fits', ignore_missing_end=True)
bad = pyfits.getdata('/home/cdeen/Data/CIAO/UT4/BadDatabase/UT4/CDMS/Common/RecnOptimiser.ModalBasis.fits', ignore_missing_end=True)


difference = good - bad

ax.imshow(difference)

fig.show()
