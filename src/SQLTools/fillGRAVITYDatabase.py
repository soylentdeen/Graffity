import sqlite3
import sys

sys.path.append('../')

import Graffity
import glob
import os

connection = sqlite3.connect(os.environ.get('GRAVITY_SQL')+'GravityObs.db')

cursor = connection.cursor()

datadir = os.environ.get('GRAVITY_DATA')

for dp, dn, fn in glob.os.walk(datadir):
    if (len(dn) == 0):
        for f in fn:
            if 'dualscip2vmred.fits' in f:
                fileBase = dp[len(datadir):]+'/'+f[:-20]
                print fileBase
                #try:
                GL = Graffity.GRAVITY_Data(fileBase=fileBase, sqlCursor=cursor)
                #GL.computeACQCAMStrehl()
                GL.addToDatabase()
                #except:
                #    print "oops!"
                #    pass



connection.commit()
connection.close()


