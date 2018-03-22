import scipy
import numpy
import astropy.io.fits as pyfits
import sys
sys.path.append('../')

import Graffity
import PlotTools
import xlrd


PlotTools.clearAllPlots()

figs, axes = PlotTools.configurePlots(1)

sheet = xlrd.open_workbook(filename='/home/grav/cdeen/GRAVITY/FibOffset_FluxCal.xlsx').sheet_by_name('S2_Flux (ABS)')

xl_names = numpy.array([sheet.col(0)[i+1].value.replace('20180224','20180309')[:-16] for i in range(sheet.nrows-1)])
xl_mjd = numpy.array([sheet.col(1)[i+1].value for i in range(sheet.nrows-1)])




print asdf
