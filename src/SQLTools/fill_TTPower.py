import sqlite3
import sys

sys.path.append('../')
import numpy
import Graffity
import glob
import os
import astropy.io.fits as pyfits

from matplotlib import pyplot

fig = pyplot.figure(0)
fig.clear()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
ax.clear()

pyplot.ion()
connection = sqlite3.connect(os.environ.get('CIAO_SQL')+'Dataloggers.db')

cursor = connection.cursor()

CDMS_BaseDir = os.environ.get('CDMS_BASEDIR')
CDMS_ConfigDir = os.environ.get('CDMS_CONFIGDIR')


command = "SELECT TIMESTAMP, DIRECTORY, TIPPOWER, TILTPOWER from CIAO_1_DataLoggers "
cursor.execute(command)
columns = cursor.fetchall()

CIAO1 = pyfits.getdata("CIAO1.fits")
CIAO2 = pyfits.getdata("CIAO2.fits")
CIAO3 = pyfits.getdata("CIAO3.fits")
CIAO4 = pyfits.getdata("CIAO4.fits")

for datalogger in columns:
    if datalogger[2] == None:
        DL = Graffity.DataLogger(directory = datalogger[1], 
                                 CDMS_BaseDir=CDMS_BaseDir,
                                 CDMS_ConfigDir=CDMS_ConfigDir, sqlCursor=cursor)
        DL.loadData()
        DL.computeStrehl(saveData=False)
        DL.computeTTPS(freq=2.2, saveData=True, Tip = CIAO1[1], Tilt = CIAO1[2])
        del(DL)

command = "SELECT TIMESTAMP, DIRECTORY, TIPPOWER, TILTPOWER from CIAO_2_DataLoggers "
cursor.execute(command)
columns = cursor.fetchall()

for datalogger in columns:
    if datalogger[2] == None:
        DL = Graffity.DataLogger(directory = datalogger[1], 
                                 CDMS_BaseDir=CDMS_BaseDir,
                                 CDMS_ConfigDir=CDMS_ConfigDir, sqlCursor=cursor)
        DL.loadData()
        DL.computeStrehl(saveData=False)
        DL.computeTTPS(freq=2.2, saveData=True, Tip = CIAO2[1], Tilt = CIAO2[2])
        del(DL)

command = "SELECT TIMESTAMP, DIRECTORY, TIPPOWER, TILTPOWER from CIAO_3_DataLoggers "
cursor.execute(command)
columns = cursor.fetchall()

for datalogger in columns:
    if datalogger[2] == None:
        DL = Graffity.DataLogger(directory = datalogger[1], 
                                 CDMS_BaseDir=CDMS_BaseDir,
                                 CDMS_ConfigDir=CDMS_ConfigDir, sqlCursor=cursor)
        DL.loadData()
        DL.computeStrehl(saveData=False)
        DL.computeTTPS(freq=2.2, saveData=True, Tip = CIAO3[1], Tilt = CIAO3[2])
        del(DL)

command = "SELECT TIMESTAMP, DIRECTORY, TIPPOWER, TILTPOWER from CIAO_4_DataLoggers "
cursor.execute(command)
columns = cursor.fetchall()

for datalogger in columns:
    if datalogger[2] == None:
        DL = Graffity.DataLogger(directory = datalogger[1], 
                                 CDMS_BaseDir=CDMS_BaseDir,
                                 CDMS_ConfigDir=CDMS_ConfigDir, sqlCursor=cursor)
        DL.loadData()
        DL.computeStrehl(saveData=False)
        DL.computeTTPS(freq=2.2, saveData=True, Tip = CIAO4[1], Tilt = CIAO4[2])
        del(DL)


connection.commit()
connection.close()


