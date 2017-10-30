import scipy
import numpy
import sys
import glob
sys.path.append('../')
import Graffity
from matplotlib import pyplot

datadir = '/home/grav/cdeen/GRAVITY/reduced/'
df = glob.glob(datadir+'GRAVI*2017-08-09T02:41:37.564_dualscip2vmred.fits')

subImageSize = 250

Telescopes= [1, 2, 3, 4]

for f in df:
    AcqCamImages = pyfits.getdata(f, ext=16)
    AcqCamData = pyfits.getdata(f, ext=17)
    AcqCamHeader = pyfits.getheader(f, ext=17)
    for im in AcqCamImages:
        for i in Telecopes:
            j = 1
            while True:
                x = AcqCamHeader.get(
