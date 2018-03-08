import sys
sys.path.append('../')
import Graffity
import numpy

filebase = '/gvstore1/forFrank/2017-07-09/reduced_20180307/GRAVI.2017-07-10T06:06:14.135'

Grav = Graffity.GRAVITY_Dual_Sci_P2VM(fileBase=filebase, processAcqCamData=True)
