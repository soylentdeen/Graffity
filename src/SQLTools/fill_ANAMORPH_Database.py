import sqlite3
import Graffity
import glob
import os
from matplotlib import pyplot

fig = pyplot.figure(0)
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

connection = sqlite3.connect(os.environ.get('CIAO_SQL')+'Dataloggers.db')

cursor = connection.cursor()

datadir = os.environ.get('CIAO_DATA')

for dp, dn, fn in glob.os.walk(datadir):
    if (len(dn) == 0) and ('ANAMORPH' in dp) and (len(glob.glob(dp+'/CIAO_LOOP_0001.fits')) > 0):
        print dp[len(datadir):]
        ANAM = Graffity.Anamorphose(directory=dp[len(datadir):], sqlCursor=cursor)
        ANAM.loadData(ax=ax)
        #DL.addToDatabase()
        print asdf



connection.commit()
connection.close()


