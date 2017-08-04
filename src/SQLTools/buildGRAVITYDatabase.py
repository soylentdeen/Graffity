import scipy
import sqlite3
import os

connection = sqlite3.connect(os.environ.get('GRAVITY_SQL')+'GravityObs.db')

cursor = connection.cursor()

cursor.execute("""DROP TABLE GRAVITY_OBS;""")

sqlCommand = """
CREATE TABLE GRAVITY_OBS ( 
TIMESTAMP FLOAT PRIMARY KEY, 
DIRECTORY VARCHAR(100),
ACQCAM_1_STREHL FLOAT,
ACQCAM_2_STREHL FLOAT, 
ACQCAM_3_STREHL FLOAT, 
ACQCAM_4_STREHL FLOAT);"""

cursor.execute(sqlCommand)

connection.commit()

connection.close()


