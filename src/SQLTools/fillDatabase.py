import sqlite3
import Graffity
import glob
import os

connection = sqlite3.connect(os.environ.get('CIAO_SQL')+'Dataloggers.db')

cursor = connection.cursor()

datadir = os.environ.get('CIAO_DATA')
CDMS_BaseDir = os.environ.get('CDMS_BASEDIR')#'/data/cdeen/CIAO_Commissioning/spcicfg'
CDMS_ConfigDir = os.environ.get('CDMS_CONFIGDIR')#'/config/RTCDATA/CIAO/DEFAULT/'

for dp, dn, fn in glob.os.walk(datadir):
    if (len(dn) == 0) and ('DATA_LOGGER' in dp) and (len(glob.glob(dp+'/CIAO_LOOP_0001.fits')) > 0):
        print dp[len(datadir):]
        DL = Graffity.DataLogger(directory=dp[len(datadir):], CDMS_BaseDir=CDMS_BaseDir,
                                 CDMS_ConfigDir=CDMS_ConfigDir,
                                 sqlCursor=cursor)
        DL.loadData()
        DL.addToDatabase()



connection.commit()
connection.close()


