import scipy
import numpy
import astropy.io.fits as pyfits
import Graffity
import glob
import matplotlib.pyplot as pyplot

datadir = '/home/cdeen/Data/CIAO/UT3/TemplateData/2016-09-01_3/FOV-212316/'

datafiles = glob.glob(datadir+'CIAO_LOOP*')

fig = pyplot.figure(0)
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

for df in datafiles:
    data = Graffity.CircularBuffer(df)
    data.processFOV(ax)



