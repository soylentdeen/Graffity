import sqlite3
import Graffity
import glob
import os

connection = sqlite3.connect('/data/cdeen/Data/CIAO/SQL/Dataloggers.db')

cursor = connection.cursor()

datadir = '/data/cdeen/Data/CIAO/CIAO_1/'
CDMS_BaseDir = '/data/cdeen/CIAO_Commissioning/spcicfg'
CDMS_ConfigDir = '/config/RTCDATA/CIAO/DEFAULT/'

for dp, dn, fn in glob.os.walk(datadir):
    if (len(dn) == 0) and ('DATA_LOGGER' in dp) and (len(glob.glob(dp+'/CIAO_LOOP_0001.fits')) > 0):
        print dp
        DL = Graffity.DataLogger(directory=dp, CDMS_BaseDir=CDMS_BaseDir, CDMS_ConfigDir=CDMS_ConfigDir,
                                 sqlCursor=cursor)
        DL.loadData()
        DL.addToDatabase()



connection.commit()
connection.close()


