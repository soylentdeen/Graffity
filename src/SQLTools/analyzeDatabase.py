import sqlite3
from matplotlib import pyplot
import scipy
import numpy
import Graffity
import os
from matplotlib import cm

CDMS_BaseDir = '/home/cdeen/Code/CIAO/SPARTA/SPARTA_CIAO/CONFIG/spcicfg'
CDMS_ConfigDir = '/config/RTCDATA/CIAO/DEFAULT/'
fig = pyplot.figure(0)
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
ax.clear()

database = '/home/cdeen/Data/CIAO/Dataloggers.db'

connection = sqlite3.connect(database)

cursor = connection.cursor()

#Get Column Names
command = "SELECT TIMESTAMP, FILENAME, ALT, AZ, PMTIP_ENC, PMTIL_ENC, STREHL, SEEING from CIAO_4_DataLoggers "
cursor.execute(command)
columns1 = cursor.fetchall()

inner_shift = []
outer_shift = []
time = []
alt = []
az = []
tip = []
tilt = []
for datalogger in columns1:
    dirname = os.path.dirname(datalogger[1])
    dl = Graffity.DataLogger(directory = dirname, CDMS_BaseDir=CDMS_BaseDir,
                             CDMS_ConfigDir=CDMS_ConfigDir)
    dl.loadData()
    dl.calculatePhotometricPupilOffset()
    inner_shift.append(dl.innerRing)
    outer_shift.append(dl.outerRing)
    alt.append(datalogger[2])
    az.append(datalogger[3])
    tip.append(datalogger[4])
    tilt.append(datalogger[5])
    time.append(datalogger[0])

time = numpy.array(time)
inner_shift= numpy.array(inner_shift)
outer_shift= numpy.array(outer_shift)
alt = numpy.array(alt)
az = numpy.array(az)
tip = numpy.array(tip)
tilt = numpy.array(tilt)

#color = (time-numpy.min(time))
#color = color/(numpy.max(time) - numpy.min(time))*128.0

r = (inner_shift[:,0]**2.0 + inner_shift[:,1]**2.0)**0.5
color = (r - numpy.min(r))
color = color/(numpy.max(r) - numpy.min(r))*256.0
#r = (outer_shift[:,0]**2.0 + outer_shift[:,1]**2.0)**0.5

#ax.scatter(time, time, c=color, cmap = 'rainbow')
#ax.scatter(tip, tilt, c=color, cmap = 'rainbow')
#ax.scatter(inner_shift[:,0], inner_shift[:,1])
#ax.scatter(outer_shift[:,0], outer_shift[:,1])
ax.scatter(alt, az, c=color, cmap='rainbow')
#ax.scatter(alt, r, c=color, cmap='rainbow')
#ax.scatter(az, r, c=color, cmap='rainbow')
#for inner, outer, t, z in zip(inner_shift, outer_shift, alt, az):
#    ax.plot([z, z], [inner[0], inner[1]], color = 'k')
#    ax.plot([t, t], [inner[0], inner[1]], color = 'r')
#    #ax.plot([inner[0], outer[0]], [inner[1], outer[1]], color = 'k')

fig.show()

cursor.close()
connection.close()
