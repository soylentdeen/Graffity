import sqlite3
import sys

sys.path.append('../')

import Graffity
import glob
import os

connection = sqlite3.connect(os.environ.get('CIAO_SQL')+'Dataloggers.db')

cursor = connection.cursor()

CDMS_BaseDir = os.environ.get('CDMS_BASEDIR')
CDMS_ConfigDir = os.environ.get('CDMS_CONFIGDIR')


command = "SELECT TIMESTAMP, DIRECTORY, TIP_RESIDUALS, TILT_RESIDUALS from CIAO_1_DataLoggers "
cursor.execute(command)
columns = cursor.fetchall()

for datalogger in columns:
    if datalogger[2] == None:
        DL = Graffity.DataLogger(directory = datalogger[1], 
                                 CDMS_BaseDir=CDMS_BaseDir,
                                 CDMS_ConfigDir=CDMS_ConfigDir, sqlCursor=cursor)
        DL.loadData()
        DL.computeTTResiduals(saveData=True)
        print datalogger[1]
        del(DL)


command = "SELECT TIMESTAMP, DIRECTORY, TIP_RESIDUALS, TILT_RESIDUALS from CIAO_2_DataLoggers "
cursor.execute(command)
columns = cursor.fetchall()

for datalogger in columns:
    if datalogger[2] == None:
        DL = Graffity.DataLogger(directory = datalogger[1], 
                                 CDMS_BaseDir=CDMS_BaseDir,
                                 CDMS_ConfigDir=CDMS_ConfigDir, sqlCursor=cursor)
        DL.loadData()
        DL.computeTTResiduals(saveData=True)
        print datalogger[1]
        del(DL)


command = "SELECT TIMESTAMP, DIRECTORY, TIP_RESIDUALS, TILT_RESIDUALS from CIAO_3_DataLoggers "
cursor.execute(command)
columns = cursor.fetchall()

for datalogger in columns:
    if datalogger[2] == None:
        DL = Graffity.DataLogger(directory = datalogger[1], 
                                 CDMS_BaseDir=CDMS_BaseDir,
                                 CDMS_ConfigDir=CDMS_ConfigDir, sqlCursor=cursor)
        DL.loadData()
        DL.computeTTResiduals(saveData=True)
        print datalogger[1]
        del(DL)


command = "SELECT TIMESTAMP, DIRECTORY, TIP_RESIDUALS, TILT_RESIDUALS from CIAO_4_DataLoggers "
cursor.execute(command)
columns = cursor.fetchall()

for datalogger in columns:
    if datalogger[2] == None:
        DL = Graffity.DataLogger(directory = datalogger[1], 
                                 CDMS_BaseDir=CDMS_BaseDir,
                                 CDMS_ConfigDir=CDMS_ConfigDir, sqlCursor=cursor)
        DL.loadData()
        DL.computeTTResiduals(saveData=True)
        print datalogger[1]
        del(DL)



connection.commit()
connection.close()


