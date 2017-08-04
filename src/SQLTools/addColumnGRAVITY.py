import sqlite3
import sys
sys.path.append('../')
import Graffity
import glob
import os

connection = sqlite3.connect(os.environ.get('GRAVITY_SQL')+'GravityObs.db')

cursor = connection.cursor()

sqlCommand = "ALTER TABLE GRAVITY_OBS DROP COLUMN ACQCAM_1_SEEING"
cursor.execute(sqlCommand)
sqlCommand = "ALTER TABLE GRAVITY_OBS DROP COLUMN ACQCAM_2_SEEING"
cursor.execute(sqlCommand)
sqlCommand = "ALTER TABLE GRAVITY_OBS DROP COLUMN ACQCAM_3_SEEING"
cursor.execute(sqlCommand)
sqlCommand = "ALTER TABLE GRAVITY_OBS DROP COLUMN ACQCAM_4_SEEING"
cursor.execute(sqlCommand)

connection.commit()
connection.close()


