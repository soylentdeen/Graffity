import sqlite3
import Graffity
import glob
import os

connection = sqlite3.connect('/home/cdeen/Data/CIAO/Dataloggers.db')

cursor = connection.cursor()

CDMS_BaseDir = '/home/cdeen/Code/CIAO/SPARTA/SPARTA_CIAO/CONFIG/spcicfg'
CDMS_ConfigDir = '/config/RTCDATA/CIAO/DEFAULT/'

sqlCommand = "SELECT * FROM CIAO_1_DataLoggers"

cursor.execute(sqlCommand)

result = cursor.fetchall()

for source in result:
    data = Graffity.DataLogger(directory=source[1], CMDS_BaseDir=CDMS_BaseDir,
                               CDMS_ConfigDir=CDMS_ConfigDir, sqlCursor=cursor)
    data.loadData()

"""
        DL = Graffity.DataLogger(directory=dp, CDMS_BaseDir=CDMS_BaseDir, CDMS_ConfigDir=CDMS_ConfigDir,
                                 sqlCursor=cursor)
        DL.loadData()
        DL.addToDatabase()
"""


#connection.commit()
connection.close()


