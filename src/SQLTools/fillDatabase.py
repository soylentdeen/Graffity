import sqlite3
import sys

sys.path.append('../')

import Graffity
import glob
import os

connection = sqlite3.connect(os.environ.get('CIAO_SQL')+'Dataloggers.db')

cursor = connection.cursor()

datadir = os.environ.get('CIAO_DATA')
CDMS_BaseDir = os.environ.get('CDMS_BASEDIR')
CDMS_ConfigDir = os.environ.get('CDMS_CONFIGDIR')

for dp, dn, fn in glob.os.walk(datadir):
    if (len(dn) == 0) and (('DATA_LOGGER' in dp) or ('DATA_EXPO' in dp)) and (len(glob.glob(dp+'/CIAO_LOOP_0001.fits')) > 0):
        print dp[len(datadir):]
        try:
            DL = Graffity.DataLogger(directory=dp[len(datadir):], CDMS_BaseDir=CDMS_BaseDir,
                                 CDMS_ConfigDir=CDMS_ConfigDir,
                                 sqlCursor=cursor)
            DL.loadData()
            DL.addToDatabase()
        except:
            pass



connection.commit()
connection.close()


