import sqlite3
import Graffity
import glob
import os

connection = sqlite3.connect(os.environ.get('CIAO_SQL')+'Dataloggers.db')

cursor = connection.cursor()

CDMS_BaseDir = os.environ.get('CDMS_BASEDIR')
CDMS_ConfigDir = os.environ.get('CDMS_CONFIGDIR')


command = "SELECT TIMESTAMP, DIRECTORY, ALT, AZ, PMTIP_ENC, PMTIL_ENC, STREHL, SEEING from CIAO_1_DataLoggers "
cursor.execute(command)
columns = cursor.fetchall()

for datalogger in columns:
    if datalogger[6] == None:
        DL = Graffity.DataLogger(directory = datalogger[1], 
                                 CDMS_BaseDir=CDMS_BaseDir,
                                 CDMS_ConfigDir=CDMS_ConfigDir, sqlCursor=cursor)
        DL.loadData()
        DL.computeStrehl(saveData=True, skipLongRecords=True)
        print datalogger[1], DL.Seeing*DL.Arcsec
        del(DL)


command = "SELECT TIMESTAMP, DIRECTORY, ALT, AZ, PMTIP_ENC, PMTIL_ENC, STREHL, SEEING from CIAO_2_DataLoggers "
cursor.execute(command)
columns = cursor.fetchall()

for datalogger in columns:
    if datalogger[6] == None:
        DL = Graffity.DataLogger(directory = datalogger[1], 
                                 CDMS_BaseDir=CDMS_BaseDir,
                                 CDMS_ConfigDir=CDMS_ConfigDir, sqlCursor=cursor)
        DL.loadData()
        DL.computeStrehl(saveData=True, skipLongRecords=True)
        print datalogger[1], DL.Seeing*DL.Arcsec
        del(DL)


command = "SELECT TIMESTAMP, DIRECTORY, ALT, AZ, PMTIP_ENC, PMTIL_ENC, STREHL, SEEING from CIAO_3_DataLoggers "
cursor.execute(command)
columns = cursor.fetchall()

for datalogger in columns:
    if datalogger[6] == None:
        DL = Graffity.DataLogger(directory = datalogger[1], 
                                 CDMS_BaseDir=CDMS_BaseDir,
                                 CDMS_ConfigDir=CDMS_ConfigDir, sqlCursor=cursor)
        DL.loadData()
        DL.computeStrehl(saveData=True, skipLongRecords=True)
        print datalogger[1], DL.Seeing*DL.Arcsec
        del(DL)


command = "SELECT TIMESTAMP, DIRECTORY, ALT, AZ, PMTIP_ENC, PMTIL_ENC, STREHL, SEEING from CIAO_4_DataLoggers "
cursor.execute(command)
columns = cursor.fetchall()

for datalogger in columns:
    if datalogger[6] == None:
        DL = Graffity.DataLogger(directory = datalogger[1], 
                                 CDMS_BaseDir=CDMS_BaseDir,
                                 CDMS_ConfigDir=CDMS_ConfigDir, sqlCursor=cursor)
        DL.loadData()
        DL.computeStrehl(saveData=True, skipLongRecords=True)
        print datalogger[1], DL.Seeing*DL.Arcsec
        del(DL)



connection.commit()
connection.close()


