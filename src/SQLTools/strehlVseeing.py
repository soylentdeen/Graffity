import sqlite3
from matplotlib import pyplot
import scipy
import numpy
import Graffity
import os
from matplotlib import cm

CDMS_BaseDir = os.environ.get('CDMS_BASEDIR')
CDMS_ConfigDir = os.environ.get('CDMS_CONFIGDIR')

fig = pyplot.figure(0)
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
ax.clear()

database = os.environ.get('CIAO_SQL')+'Dataloggers.db'

connection = sqlite3.connect(database)

cursor = connection.cursor()

Strehl = {}
Seeing = {}

colors = ['b', 'g', 'r', 'y']
for i in numpy.arange(4):
    ciaoID = i+1
    #Get Column Names
    command = "SELECT STREHL, SEEING, WINDSP from CIAO_%d_DataLoggers " % ciaoID
    
    cursor.execute(command)
    data = numpy.array(cursor.fetchall())

    Strehl[ciaoID] = data[:,0]
    Seeing[ciaoID] = data[:,1]

    ax.scatter(data[:,2], data[:,0], color = colors[i])

fig.show()

cursor.close()
connection.close()
