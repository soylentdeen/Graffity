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


"""
command = "SELECT TIMESTAMP, DIRECTORY, ALT, AZ, PMTIP_ENC, PMTIL_ENC, STREHL, SEEING from CIAO_1_DataLoggers "
cursor.execute(command)
columns = cursor.fetchall()

Tips = []
Tilts = []

for datalogger in columns:
    if datalogger[0] != None:
        DL = Graffity.DataLogger(directory = datalogger[1], 
                                 CDMS_BaseDir=CDMS_BaseDir,
                                 CDMS_ConfigDir=CDMS_ConfigDir, sqlCursor=cursor)
        DL.loadData()
        DL.computeStrehl(saveData=False)
        Tip, Tilt, Freq = DL.computeTTPS(freq=2.2, saveData=False, returnTipTilt=True)
        Tips.append(Tip)
        Tilts.append(Tilt)
        del(DL)


Tips = numpy.median(numpy.array(Tips), axis=0)
Tilts = numpy.median(numpy.array(Tilts), axis=0)

HDU = pyfits.PrimaryHDU(data=[Freq, Tips, Tilts])
HDU.writeto("CIAO1.fits")

"""

command = "SELECT TIMESTAMP, DIRECTORY, ALT, AZ, PMTIP_ENC, PMTIL_ENC, STREHL, SEEING from CIAO_2_DataLoggers "
cursor.execute(command)
columns = cursor.fetchall()

Tips = []
Tilts = []

for datalogger in columns:
    if datalogger[0] != None:
        DL = Graffity.DataLogger(directory = datalogger[1], 
                                 CDMS_BaseDir=CDMS_BaseDir,
                                 CDMS_ConfigDir=CDMS_ConfigDir, sqlCursor=cursor)
        DL.loadData()
        DL.computeStrehl(saveData=False)
        Tip, Tilt, Freq = DL.computeTTPS(freq=2.2, saveData=False, returnTipTilt=True)
        Tips.append(Tip)
        Tilts.append(Tilt)
        del(DL)

Tips = numpy.median(numpy.array(Tips), axis=0)
Tilts = numpy.median(numpy.array(Tilts), axis=0)

HDU = pyfits.PrimaryHDU(data=[Freq, Tips, Tilts])
HDU.writeto("CIAO2.fits")

command = "SELECT TIMESTAMP, DIRECTORY, ALT, AZ, PMTIP_ENC, PMTIL_ENC, STREHL, SEEING from CIAO_3_DataLoggers "
cursor.execute(command)
columns = cursor.fetchall()

Tips = []
Tilts = []

for datalogger in columns:
    if datalogger[0] != None:
        DL = Graffity.DataLogger(directory = datalogger[1], 
                                 CDMS_BaseDir=CDMS_BaseDir,
                                 CDMS_ConfigDir=CDMS_ConfigDir, sqlCursor=cursor)
        DL.loadData()
        DL.computeStrehl(saveData=False)
        Tip, Tilt, Freq = DL.computeTTPS(freq=2.2, saveData=False, returnTipTilt=True)
        Tips.append(Tip)
        Tilts.append(Tilt)
        del(DL)

Tips = numpy.median(numpy.array(Tips), axis=0)
Tilts = numpy.median(numpy.array(Tilts), axis=0)

HDU = pyfits.PrimaryHDU(data=[Freq, Tips, Tilts])
HDU.writeto("CIAO3.fits")

command = "SELECT TIMESTAMP, DIRECTORY, ALT, AZ, PMTIP_ENC, PMTIL_ENC, STREHL, SEEING from CIAO_4_DataLoggers "
cursor.execute(command)
columns = cursor.fetchall()

Tips = []
Tilts = []

for datalogger in columns:
    if datalogger[0] != None:
        DL = Graffity.DataLogger(directory = datalogger[1], 
                                 CDMS_BaseDir=CDMS_BaseDir,
                                 CDMS_ConfigDir=CDMS_ConfigDir, sqlCursor=cursor)
        DL.loadData()
        DL.computeStrehl(saveData=False)
        Tip, Tilt, Freq = DL.computeTTPS(freq=2.2, saveData=False, returnTipTilt=True)
        Tips.append(Tip)
        Tilts.append(Tilt)
        del(DL)


Tips = numpy.median(numpy.array(Tips), axis=0)
Tilts = numpy.median(numpy.array(Tilts), axis=0)

HDU = pyfits.PrimaryHDU(data=[Freq, Tips, Tilts])
HDU.writeto("CIAO4.fits")

connection.commit()
connection.close()


