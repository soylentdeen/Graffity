import scipy
import numpy
import Graffity
import matplotlib.pyplot as pyplot
import sys
import glob

fig = pyplot.figure(0)
fig.clear()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])


datadir = '/home/deen/Data/GRAVITY/MATLAB_DATA/'

date = '2016-03-29/'
template = sys.argv[1]

loopFiles = glob.glob(datadir+date+template+'/CIAO_LOOP_*.fits')

for loop in loopFiles:
    CB = Graffity.CircularBuffer(loop)
    CB.loadTemplateMaps()

    CB.processSlopes()
    CB.processVoltages()

    #ax.bar(numpy.arange(len(CB.ProcessedVoltages.zerns_average)),
    #        CB.ProcessedVoltages.zerns_average, width=0.4)
    #ax.bar(numpy.arange(len(CB.ProcessedVoltages.zerns_std))+0.5,
    #        CB.ProcessedVoltages.zerns_std, width=0.4)
    #ax.plot(CB.ProcessedGradients.zerns)
    ax.plot(CB.ProcessedGradients.average)

    # difference to DM Flat
    #   - per mode
    #   - per Zernike
    # Saturations
    #   - % per frame
    #   - rate per actuator

    
fig.show()
