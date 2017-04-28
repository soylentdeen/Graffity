import sqlite3
from matplotlib import pyplot
import scipy
import numpy
import Graffity
import os
from matplotlib import cm

CDMS_BaseDir = os.environ.get('CDMS_BASEDIR')
CDMS_ConfigDir = os.environ.get('CDMS_CONFIGDIR')

fig1 = pyplot.figure(0)
ax11 = fig1.add_axes([0.1, 0.1, 0.4, 0.4])
ax12 = fig1.add_axes([0.1, 0.4, 0.4, 0.4])
ax13 = fig1.add_axes([0.4, 0.1, 0.4, 0.4])
ax14 = fig1.add_axes([0.4, 0.4, 0.4, 0.4])
ax11.clear()
ax12.clear()
ax13.clear()
ax14.clear()

database = os.environ.get('CIAO_SQL')+'Dataloggers.db'

connection = sqlite3.connect(database)

cursor = connection.cursor()

Strehl = {}
Seeing = {}

colors = ['b', 'g', 'r', 'y']
for i in numpy.arange(4):
    ciaoID = i+1
    #Get Column Names
    command = "SELECT TIMESTAMP, FLDLX, FLDLY, FSM_A_U, FSM_A_W, FSM_A_X, FSM_A_Y, FSM_B_U, FSM_B_W, FSM_B_X, FSM_B_Y, VCM_A_U, VCM_A_W, VCM_A_X, VCM_A_Y, VCM_B_U, VCM_B_W, VCM_B_X, VCM_B_Y, M10_POSANG from CIAO_%d_DataLoggers " % ciaoID
    
    cursor.execute(command)
    data = numpy.array(cursor.fetchall())

    ax11.scatter(data[:,1], data[:,2], color = colors[i])
    ax12.scatter(data[:,3], data[:,4], color = colors[i])

fig1.show()

cursor.close()
connection.close()
