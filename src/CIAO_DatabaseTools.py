import sqlite3
import numpy
import time
from datetime import datetime
import ephem
import os
from astropy import time as aptime

class GRAVITY_Database( object ):
    def __init__(self, name = 'GRAVITY.db'):
        self.database = os.environ.get('GRAVITY_SQL')+name
        self.connection = sqlite3.connect(self.database)
        self.cursor = self.connection.cursor()
        self.Paranal = ephem.Observer()
        self.Paranal.lat=-24.62694
        self.Paranal.long=-70.405
        self.Sun = ephem.Sun(self.Paranal)

    def addTable(self, sqlCommand='', tableName='GRAVITY_REDUCED'):
        if sqlCommand == '':
            sqlCommand = """
            CREATE TABLE %s ( 
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
            HOCTR_PRA_GAIN FLOAT);""" % tableName
        self.cursor.execute(sqlCommand)

    def query(self, keywords={}, timeOfDay='BOTH', startTime=None, endTime=None):
        retval = []
        sqlCommand = "SELECT"
        first = True
        for keyword in keywords:
            if first:
                first = False
            else:
                sqlCommand = sqlCommand + ", "
            sqlCommand = sqlCommand + keyword

        sqlCommand = sqlCommand + ", TIMESTAMP, "
        sqlCommand = sqlCommand + " from GRAVITY_REDUCED "

        self.cursor.execute(sqlCommand)
        result = numpy.array(self.cursor.fetchall())

        if timeOfDay == 'NIGHT':  #Night
            nightTime = []
            for observation in result:
                if self.atNight(observation[-4]):
                    nightTime.append(observation)
            result = numpy.array(nightTime)
        elif timeOfDay == 'DAY': #Day
            dayTime = []
            for observation in result:
                if self.atNight(observation[-4]):
                    dayTime.append(observation)
            result = numpy.array(dayTime)

        if startTime != None:
            result = result[result[:,-4] > time.mktime(time.strptime(startTime, '%Y-%m-%d %H:%M:%S'))]
        if endTime != None:
            result = result[result[:,-4] < time.mktime(time.strptime(endTime, '%Y-%m-%d %H:%M:%S'))]
        retval[UT] = result

        return retval

    def atNight(self, timestamp):
        ts = datetime.utcfromtimestamp(float(timestamp))
        previous_rising = ephem.Date(self.Paranal.previous_rising(self.Sun, start=ts))
        previous_setting = ephem.Date(self.Paranal.previous_setting(self.Sun, start=ts))
        current_time = ephem.Date(ts)
        
        if (current_time - previous_rising) > (current_time - previous_setting):
            return True
        else:
            return False
        
    def close(self):
        self.connection.close()

class CIAO_Database( object ):
    def __init__(self, name='Dataloggers.db'):
        self.database = os.environ.get('CIAO_SQL')+name
        self.connection = sqlite3.connect(self.database)
        self.cursor = self.connection.cursor()
        self.Paranal = ephem.Observer()
        self.Paranal.lat=-24.62694
        self.Paranal.long=-70.405
        self.Sun = ephem.Sun(self.Paranal)
    
    def addTable(self, columnNames = "", tableName=''):
        sqlCommand = "CREATE TABLE IF NOT EXISTS %s (%s);" % (tableName, columnNames)
        self.cursor.execute(sqlCommand)

    def query(self, keywords = {}, UTS=[1, 2, 3, 4], timeOfDay='BOTH', startTime=None, endTime=None,
              AVC_State='EITHER', Pupil_State=0, TemplateType='DataLoggers'):

        retval = {}
        for UT in UTS:
            sqlCommand = "SELECT "
            first = True
            for keyword in keywords:
                if first:
                    first = False
                else:
                    sqlCommand = sqlCommand + ", "
                sqlCommand = sqlCommand + keyword
        
            sqlCommand = sqlCommand + ", TIMESTAMP, DIRECTORY, AVC_STATE, PUP_TRK" 
            sqlCommand = sqlCommand + " from CIAO_"+str(UT)+"_%s " %TemplateType
        
            self.cursor.execute(sqlCommand)
            result = numpy.array(self.cursor.fetchall())
    
            if AVC_State == 'ON':
                result = result[result[:,-2] == u'True']
            elif AVC_State == 'OFF':
                result = result[result[:,-2] == u'False']
    
            if timeOfDay == 'NIGHT':  #Night
                nightTime = []
                for observation in result:
                    if self.atNight(observation[-4]):
                        nightTime.append(observation)
                result = numpy.array(nightTime)
            elif timeOfDay == 'DAY':  #Day
                dayTime = []
                for observation in result:
                    if not(self.atNight(observation[-4])):
                        dayTime.append(observation)
                result = numpy.array(dayTime)

            if startTime != None:
                result = result[result[:,-4] > aptime.Time(startTime, format='iso').mjd]
                #result = result[result[:,-4] > time.mktime(time.strptime(startTime, '%Y-%m-%d %H:%M:%S'))]
            if endTime != None:
                result = result[result[:,-4] < aptime.Time(endTime, format='iso').mjd]
                #result = result[result[:,-4] < time.mktime(time.strptime(endTime, '%Y-%m-%d %H:%M:%S'))]
            retval[UT] = result
        return retval

    def atNight(self, timestamp):
        #ts = datetime.utcfromtimestamp(float(timestamp))
        ts = aptime.Time(timestamp, format='mjd')   #This is obviously wrong
        previous_rising = ephem.Date(self.Paranal.previous_rising(self.Sun, start=ts.to_datetime()))
        previous_setting = ephem.Date(self.Paranal.previous_setting(self.Sun, start=ts.to_datetime()))
        current_time = ephem.Date(ts.to_datetime())
        
        if (current_time - previous_rising) > (current_time - previous_setting):
            return True
        else:
            return False
        
    def close(self):
        self.connection.close()
