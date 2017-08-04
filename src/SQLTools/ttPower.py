import sys

sys.path.append('../')

import CIAO_DatabaseTools
from matplotlib import pyplot
import numpy
import scipy


fig = pyplot.figure(0)
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
ax.clear()


db = CIAO_DatabaseTools.CIAO_Database()

UTS = [1,2,3,4]
colors = ['b', 'g', 'r', 'y']

values = db.query(keywords = ["TIP_RESIDUALS", "TILT_RESIDUALS", "TIPPOWER", "TILTPOWER"], AVC_State='ON',
        timeOfDay='NIGHT', UTS=UTS, startTime='2017-07-10 00:00:00')
                  #startTime='2017-03-25 00:00:00')
                  #endTime='2017-02-28 00:00:00')

for i in UTS:
    ax.scatter(values[i][:,4], values[i][:,2], color = colors[i-1])
    ax.scatter(values[i][:,4], values[i][:,3], color = colors[i-1])
    #ax.scatter(values[i][:,2], values[i][:,3], color = colors[i-1])

fig.show()


db.close()
