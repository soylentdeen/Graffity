import scipy
import sqlite3
import os

connection = sqlite3.connect(os.environ.get('CIAO_SQL')+'Dataloggers.db')

cursor = connection.cursor()

cursor.execute("""DROP TABLE CIAO_1_Dataloggers;""")
cursor.execute("""DROP TABLE CIAO_2_Dataloggers;""")
cursor.execute("""DROP TABLE CIAO_3_Dataloggers;""")
cursor.execute("""DROP TABLE CIAO_4_Dataloggers;""")

sqlCommand = """
CREATE TABLE CIAO_1_DataLoggers ( 
TIMESTAMP FLOAT PRIMARY KEY, 
DIRECTORY VARCHAR(100),
SEEING FLOAT,
ASM_SEEING FLOAT, 
STREHL FLOAT,
TAU0 FLOAT,
TERR FLOAT,
AVC_STATE BOOLEAN,
CM_MODES INTEGER,
WFS_GEOM CHARACTER(8),
GAIN FLOAT,
HOCTR_AWF_ENABLE BOOLEAN,
HOCTR_GARBAGE_GAIN FLOAT,
HOCTR_KI FLOAT,
HOCTR_KT FLOAT,
HOCTR_PRA_ENABLE BOOLEAN,
HOCTR_PRA_GAIN FLOAT,
HOCTR_SMA_ENABLE BOOLEAN,
HOCTR_SMA_HIGH FLOAT,
HOCTR_SMA_ITERATIONS INTEGER,
HOCTR_SMA_LOW FLOAT,
HOCTR_TT_KI FLOAT,
HOCTR_TT_KT FLOAT,
LOOPRATE FLOAT,
TT_LOOP_STATE BOOLEAN,
TTX_REFPOS FLOAT,
TTY_REFPOS FLOAT,
VIB_SR FLOAT,
WFS_MODE CHARACTER(8),
IM_TRK_MODE CHARACTER(8),
ALT FLOAT,
AZ FLOAT,
DIT FLOAT,
DROT_ENC FLOAT,
FILT_ENC FLOAT,
MSEL_ENC FLOAT,
MSEL_NAME FLOAT,
PMTIL_ENC FLOAT,
PMTIP_ENC FLOAT,
FLDLX FLOAT,
FLDLY FLOAT,
FSM_A_U FLOAT,
FSM_A_W FLOAT,
FSM_A_X FLOAT,
FSM_A_Y FLOAT,
FSM_B_U FLOAT,
FSM_B_W FLOAT,
FSM_B_X FLOAT,
FSM_B_Y FLOAT,
VCM_A_U FLOAT,
VCM_A_W FLOAT,
VCM_A_X FLOAT,
VCM_A_Y FLOAT,
VCM_B_U FLOAT,
VCM_B_W FLOAT,
VCM_B_X FLOAT,
VCM_B_Y FLOAT,
M10_POSANG FLOAT,
RHUM FLOAT,
TEMP FLOAT,
THETA0 FLOAT,
WINDDIR FLOAT,
WINDSP FLOAT,
PRLTIC FLOAT,
PUP_TRK BOOLEAN,
FLUX FLOAT,
TIP_RESIDUALS FLOAT,
TILT_RESIDUALS FLOAT);"""

cursor.execute(sqlCommand)

sqlCommand = """
CREATE TABLE CIAO_2_DataLoggers ( 
TIMESTAMP FLOAT PRIMARY KEY, 
DIRECTORY VARCHAR(100),
SEEING FLOAT,
ASM_SEEING FLOAT, 
STREHL FLOAT,
TAU0 FLOAT,
TERR FLOAT,
AVC_STATE BOOLEAN,
CM_MODES INTEGER,
WFS_GEOM CHARACTER(8),
GAIN FLOAT,
HOCTR_AWF_ENABLE BOOLEAN,
HOCTR_GARBAGE_GAIN FLOAT,
HOCTR_KI FLOAT,
HOCTR_KT FLOAT,
HOCTR_PRA_ENABLE BOOLEAN,
HOCTR_PRA_GAIN FLOAT,
HOCTR_SMA_ENABLE BOOLEAN,
HOCTR_SMA_HIGH FLOAT,
HOCTR_SMA_ITERATIONS INTEGER,
HOCTR_SMA_LOW FLOAT,
HOCTR_TT_KI FLOAT,
HOCTR_TT_KT FLOAT,
LOOPRATE FLOAT,
TT_LOOP_STATE BOOLEAN,
TTX_REFPOS FLOAT,
TTY_REFPOS FLOAT,
VIB_SR FLOAT,
WFS_MODE CHARACTER(8),
IM_TRK_MODE CHARACTER(8),
ALT FLOAT,
AZ FLOAT,
DIT FLOAT,
DROT_ENC FLOAT,
FILT_ENC FLOAT,
MSEL_ENC FLOAT,
MSEL_NAME FLOAT,
PMTIL_ENC FLOAT,
PMTIP_ENC FLOAT,
FLDLX FLOAT,
FLDLY FLOAT,
FSM_A_U FLOAT,
FSM_A_W FLOAT,
FSM_A_X FLOAT,
FSM_A_Y FLOAT,
FSM_B_U FLOAT,
FSM_B_W FLOAT,
FSM_B_X FLOAT,
FSM_B_Y FLOAT,
VCM_A_U FLOAT,
VCM_A_W FLOAT,
VCM_A_X FLOAT,
VCM_A_Y FLOAT,
VCM_B_U FLOAT,
VCM_B_W FLOAT,
VCM_B_X FLOAT,
VCM_B_Y FLOAT,
M10_POSANG FLOAT,
RHUM FLOAT,
TEMP FLOAT,
THETA0 FLOAT,
WINDDIR FLOAT,
WINDSP FLOAT,
PRLTIC FLOAT,
PUP_TRK BOOLEAN,
FLUX FLOAT,
TIP_RESIDUALS FLOAT,
TILT_RESIDUALS FLOAT);"""

cursor.execute(sqlCommand)

sqlCommand = """
CREATE TABLE CIAO_3_DataLoggers ( 
TIMESTAMP FLOAT PRIMARY KEY, 
DIRECTORY VARCHAR(100),
SEEING FLOAT,
ASM_SEEING FLOAT, 
STREHL FLOAT,
TAU0 FLOAT,
TERR FLOAT,
AVC_STATE BOOLEAN,
CM_MODES INTEGER,
WFS_GEOM CHARACTER(8),
GAIN FLOAT,
HOCTR_AWF_ENABLE BOOLEAN,
HOCTR_GARBAGE_GAIN FLOAT,
HOCTR_KI FLOAT,
HOCTR_KT FLOAT,
HOCTR_PRA_ENABLE BOOLEAN,
HOCTR_PRA_GAIN FLOAT,
HOCTR_SMA_ENABLE BOOLEAN,
HOCTR_SMA_HIGH FLOAT,
HOCTR_SMA_ITERATIONS INTEGER,
HOCTR_SMA_LOW FLOAT,
HOCTR_TT_KI FLOAT,
HOCTR_TT_KT FLOAT,
LOOPRATE FLOAT,
TT_LOOP_STATE BOOLEAN,
TTX_REFPOS FLOAT,
TTY_REFPOS FLOAT,
VIB_SR FLOAT,
WFS_MODE CHARACTER(8),
IM_TRK_MODE CHARACTER(8),
ALT FLOAT,
AZ FLOAT,
DIT FLOAT,
DROT_ENC FLOAT,
FILT_ENC FLOAT,
MSEL_ENC FLOAT,
MSEL_NAME FLOAT,
PMTIL_ENC FLOAT,
PMTIP_ENC FLOAT,
FLDLX FLOAT,
FLDLY FLOAT,
FSM_A_U FLOAT,
FSM_A_W FLOAT,
FSM_A_X FLOAT,
FSM_A_Y FLOAT,
FSM_B_U FLOAT,
FSM_B_W FLOAT,
FSM_B_X FLOAT,
FSM_B_Y FLOAT,
VCM_A_U FLOAT,
VCM_A_W FLOAT,
VCM_A_X FLOAT,
VCM_A_Y FLOAT,
VCM_B_U FLOAT,
VCM_B_W FLOAT,
VCM_B_X FLOAT,
VCM_B_Y FLOAT,
M10_POSANG FLOAT,
RHUM FLOAT,
TEMP FLOAT,
THETA0 FLOAT,
WINDDIR FLOAT,
WINDSP FLOAT,
PRLTIC FLOAT,
PUP_TRK BOOLEAN,
FLUX FLOAT,
TIP_RESIDUALS FLOAT,
TILT_RESIDUALS FLOAT);"""

cursor.execute(sqlCommand)

sqlCommand = """
CREATE TABLE CIAO_4_DataLoggers ( 
TIMESTAMP FLOAT PRIMARY KEY, 
DIRECTORY VARCHAR(100),
SEEING FLOAT,
ASM_SEEING FLOAT, 
STREHL FLOAT,
TAU0 FLOAT,
TERR FLOAT,
AVC_STATE BOOLEAN,
CM_MODES INTEGER,
WFS_GEOM CHARACTER(8),
GAIN FLOAT,
HOCTR_AWF_ENABLE BOOLEAN,
HOCTR_GARBAGE_GAIN FLOAT,
HOCTR_KI FLOAT,
HOCTR_KT FLOAT,
HOCTR_PRA_ENABLE BOOLEAN,
HOCTR_PRA_GAIN FLOAT,
HOCTR_SMA_ENABLE BOOLEAN,
HOCTR_SMA_HIGH FLOAT,
HOCTR_SMA_ITERATIONS INTEGER,
HOCTR_SMA_LOW FLOAT,
HOCTR_TT_KI FLOAT,
HOCTR_TT_KT FLOAT,
LOOPRATE FLOAT,
TT_LOOP_STATE BOOLEAN,
TTX_REFPOS FLOAT,
TTY_REFPOS FLOAT,
VIB_SR FLOAT,
WFS_MODE CHARACTER(8),
IM_TRK_MODE CHARACTER(8),
ALT FLOAT,
AZ FLOAT,
DIT FLOAT,
DROT_ENC FLOAT,
FILT_ENC FLOAT,
MSEL_ENC FLOAT,
MSEL_NAME FLOAT,
PMTIL_ENC FLOAT,
PMTIP_ENC FLOAT,
FLDLX FLOAT,
FLDLY FLOAT,
FSM_A_U FLOAT,
FSM_A_W FLOAT,
FSM_A_X FLOAT,
FSM_A_Y FLOAT,
FSM_B_U FLOAT,
FSM_B_W FLOAT,
FSM_B_X FLOAT,
FSM_B_Y FLOAT,
VCM_A_U FLOAT,
VCM_A_W FLOAT,
VCM_A_X FLOAT,
VCM_A_Y FLOAT,
VCM_B_U FLOAT,
VCM_B_W FLOAT,
VCM_B_X FLOAT,
VCM_B_Y FLOAT,
M10_POSANG FLOAT,
RHUM FLOAT,
TEMP FLOAT,
THETA0 FLOAT,
WINDDIR FLOAT,
WINDSP FLOAT,
PRLTIC FLOAT,
PUP_TRK BOOLEAN,
FLUX FLOAT,
TIP_RESIDUALS FLOAT,
TILT_RESIDUALS FLOAT);"""

cursor.execute(sqlCommand)

connection.commit()

connection.close()


