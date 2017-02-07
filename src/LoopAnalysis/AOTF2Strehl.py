import scipy
import numpy
from matplotlib import pyplot
import astropy.io.fits as pyfits
import Graffity

fig = pyplot.figure(0)
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

for i in range(4):
    df = '/media/cdeen/My Passport/GRAVITY/CIAO#1/Paranal/2016-04-03/Perf_vs_Magn/AO_TF-092142/CIAO_LOOP_000'+str(i+1)+'.fits'

    data = pyfits.getdata(df)

    slopes = data.field('Gradients')
    HODM = data.field('HODM_Positions')
    header = pyfits.getheader(df)
    print header.get('ESO AOS GLOBAL GAIN')

ax.plot(slopes[:,0])
fig.show()
