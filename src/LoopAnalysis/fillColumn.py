import sqlite3
import Graffity
import glob
import os

connection = sqlite3.connect('/home/cdeen/Data/CIAO/Dataloggers.db')

cursor = connection.cursor()

datadir = '/home/cdeen/Data/CIAO/CIAO_2/'
CDMS_BaseDir = '/home/cdeen/Code/CIAO/SPARTA/SPARTA_CIAO/CONFIG/spcicfg'
CDMS_ConfigDir = '/config/RTCDATA/CIAO/DEFAULT/'

for dp, dn, fn in glob.os.walk(datadir):
    if (len(dn) == 0) and ('DATA_LOGGER' in dp):
        print dp
        DL = Graffity.DataLogger(directory=dp, CDMS_BaseDir=CDMS_BaseDir, CDMS_ConfigDir=CDMS_ConfigDir,
                                 sqlCursor=cursor)
        DL.loadData()
        DL.addToDatabase()



connection.commit()
connection.close()


