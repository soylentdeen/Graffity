import sqlite3
import Graffity
import glob
import os

connection = sqlite3.connect(os.environ.get('CIAO_SQL')+'Dataloggers.db')

cursor = connection.cursor()

CDMS_BaseDir = '/home/cdeen/Code/CIAO/SPARTA/SPARTA_CIAO/CONFIG/spcicfg'
CDMS_ConfigDir = '/config/RTCDATA/CIAO/DEFAULT/'

sqlCommand = "ALTER TABLE CIAO_1_DataLoggers ADD TIPPOWER"
cursor.execute(sqlCommand)
sqlCommand = "ALTER TABLE CIAO_1_DataLoggers ADD TILTPOWER"
cursor.execute(sqlCommand)

sqlCommand = "ALTER TABLE CIAO_2_DataLoggers ADD TIPPOWER"
cursor.execute(sqlCommand)
sqlCommand = "ALTER TABLE CIAO_2_DataLoggers ADD TILTPOWER"
cursor.execute(sqlCommand)

sqlCommand = "ALTER TABLE CIAO_3_DataLoggers ADD TIPPOWER"
cursor.execute(sqlCommand)
sqlCommand = "ALTER TABLE CIAO_3_DataLoggers ADD TILTPOWER"
cursor.execute(sqlCommand)

sqlCommand = "ALTER TABLE CIAO_4_DataLoggers ADD TIPPOWER"
cursor.execute(sqlCommand)
sqlCommand = "ALTER TABLE CIAO_4_DataLoggers ADD TILTPOWER"
cursor.execute(sqlCommand)


connection.commit()
connection.close()


