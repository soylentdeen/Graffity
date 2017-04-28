import scipy
import Graffity
import numpy

datadir = 'CIAO_1/2017-03/2017-03-27/DATA_LOGGER-092621'

DL = Graffity.DataLogger(directory=datadir)
DL.loadData()
DL.computeStrehl()
print DL.Strehl, DL.Seeing*DL.Arcsec
