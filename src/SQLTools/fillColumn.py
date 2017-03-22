import sqlite3
import Graffity
import glob
import os

connection = sqlite3.connect('/data/cdeen/Data/CIAO/SQL/Dataloggers.db')

cursor = connection.cursor()

CDMS_BaseDir = '/data/cdeen/CIAO_Commissioning/spcicfg'
CDMS_ConfigDir = '/config/RTCDATA/CIAO/DEFAULT/'


command = "SELECT TIMESTAMP, DIRECTORY, ALT, AZ, PMTIP_ENC, PMTIL_ENC, STREHL, SEEING from CIAO_1_DataLoggers "
cursor.execute(command)
columns = cursor.fetchall()

for datalogger in columns:
    DL = Graffity.DataLogger(directory = datalogger[1], 
                             CDMS_BaseDir=CDMS_BaseDir,
                             CDMS_ConfigDir=CDMS_ConfigDir, sqlCursor=cursor)
    DL.loadData()
    DL.computeStrehl(saveData=True)
    print "CIAO 1 :", datalogger[0], DL.Seeing*DL.Arcsec
    del(DL)


command = "SELECT TIMESTAMP, DIRECTORY, ALT, AZ, PMTIP_ENC, PMTIL_ENC, STREHL, SEEING from CIAO_2_DataLoggers "
cursor.execute(command)
columns = cursor.fetchall()

for datalogger in columns:
    DL = Graffity.DataLogger(directory = datalogger[1], 
                             CDMS_BaseDir=CDMS_BaseDir,
                             CDMS_ConfigDir=CDMS_ConfigDir, sqlCursor=cursor)
    DL.loadData()
    DL.computeStrehl(saveData=True)
    print "CIAO 2 :", datalogger[0], DL.Seeing*DL.Arcsec
    del(DL)


command = "SELECT TIMESTAMP, DIRECTORY, ALT, AZ, PMTIP_ENC, PMTIL_ENC, STREHL, SEEING from CIAO_3_DataLoggers "
cursor.execute(command)
columns = cursor.fetchall()

for datalogger in columns:
    DL = Graffity.DataLogger(directory = datalogger[1], 
                             CDMS_BaseDir=CDMS_BaseDir,
                             CDMS_ConfigDir=CDMS_ConfigDir, sqlCursor=cursor)
    DL.loadData()
    DL.computeStrehl(saveData=True)
    print "CIAO 3 :", datalogger[0], DL.Seeing*DL.Arcsec
    del(DL)


command = "SELECT TIMESTAMP, DIRECTORY, ALT, AZ, PMTIP_ENC, PMTIL_ENC, STREHL, SEEING from CIAO_4_DataLoggers "
cursor.execute(command)
columns = cursor.fetchall()

for datalogger in columns:
    DL = Graffity.DataLogger(directory = datalogger[1], 
                             CDMS_BaseDir=CDMS_BaseDir,
                             CDMS_ConfigDir=CDMS_ConfigDir, sqlCursor=cursor)
    DL.loadData()
    DL.computeStrehl(saveData=True)
    print "CIAO 4 :", datalogger[0], DL.Seeing*DL.Arcsec
    del(DL)



connection.commit()
connection.close()


