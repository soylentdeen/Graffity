import sqlite3
import numpy
import time
from datetime import datetime
import ephem
import os
from astropy import time as aptime

#class QueryResult( object ):
#    def __init__(self, keywords = {}):
#        self.

class GRAVITY_Database( object ):
    def __init__(self, name = 'GravityObs.db'):
        self.database = os.environ.get('GRAVITY_SQL')+name
        self.connection = sqlite3.connect(self.database)
        self.cursor = self.connection.cursor()
        self.Paranal = ephem.Observer()
        self.Paranal.lat=-24.62694
        self.Paranal.long=-70.405
        self.Sun = ephem.Sun(self.Paranal)

    def addTable(self, sqlCommand='', tableName='GRAVITY_OBS'):
        if sqlCommand == '':
            sqlCommand = """
            CREATE TABLE %s ( 
            TIMESTAMP FLOAT PRIMARY KEY, 
            DIRECTORY VARCHAR(100),
            FTOBJ_NAME VARCHAR(20),
            FTMAG FLOAT,
            SOBJ_NAME VARCHAR(20),
            SOBJMAG FLOAT
            SOBJ_SWAP BOOLEAN);""" % tableName
        self.cursor.execute(sqlCommand)

    def query(self, keywords={}, timeOfDay='BOTH', startTime=None, endTime=None):
        retval = []
        sqlCommand = "SELECT "
        first = True
        for keyword in keywords:
            if first:
                first = False
            else:
                sqlCommand = sqlCommand + ", "
            sqlCommand = sqlCommand + keyword

        if len(keywords) == 0:
            sqlCommand = sqlCommand + " TIMESTAMP, DIRECTORY"
        else:
            sqlCommand = sqlCommand + ", TIMESTAMP, DIRECTORY"
        sqlCommand = sqlCommand + " from GRAVITY_OBS "

        self.cursor.execute(sqlCommand)
        result = numpy.array(self.cursor.fetchall())

        if timeOfDay == 'NIGHT':  #Night
            nightTime = []
            for observation in result:
                if self.atNight(observation[-2]):
                    nightTime.append(observation)
            result = numpy.array(nightTime)
        elif timeOfDay == 'DAY': #Day
            dayTime = []
            for observation in result:
                if self.atNight(observation[-2]):
                    dayTime.append(observation)
            result = numpy.array(dayTime)

        if startTime != None:
            result = result[result[:,-2] > unicode(aptime.Time(startTime, format='iso').mjd)]
        if endTime != None:
            result = result[result[:,-2] < unicode(aptime.Time(endTime, format='iso').mjd)]

        #retval = []
        #fullKeywords = keywords + ['TIMESTAMP', 'DIRECTORY']
        #for r in result:
        #    retDict = {}
        #    for key, val in zip(fullKeywords, r):
        #        retDict[key] = val
        #    retval.append(retDict)
        return result

    def atNight(self, timestamp):
        ts = aptime.Time(float(timestamp), format='mjd')
        previous_rising = ephem.Date(self.Paranal.previous_rising(self.Sun,
            start=ts.to_datetime()))
        previous_setting = ephem.Date(self.Paranal.previous_setting(self.Sun,
            start=ts.to_datetime()))
        current_time = ephem.Date(ts.to_datetime())
        
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
        fullKeywords = keywords + ['TIMESTAMP', 'DIRECTORY']
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
        
            print sqlCommand
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
            if endTime != None:
                result = result[result[:,-4] < aptime.Time(endTime, format='iso').mjd]
            #retval[UT] = []
            #for r in result:
            #    retDict = {}
            #    for key, val in zip(fullKeywords, r):
            #        retDict[key] = val
            #    retval[UT].append(retDict)
            retval[UT] = result
        return retval

        return retval

    def atNight(self, timestamp):
        #ts = datetime.utcfromtimestamp(float(timestamp))
        ts = aptime.Time(float(timestamp), format='mjd')   #This is obviously wrong
        previous_rising = ephem.Date(self.Paranal.previous_rising(self.Sun, start=ts.to_datetime()))
        previous_setting = ephem.Date(self.Paranal.previous_setting(self.Sun, start=ts.to_datetime()))
        current_time = ephem.Date(ts.to_datetime())
        
        if (current_time - previous_rising) > (current_time - previous_setting):
            return True
        else:
            return False
        
    def close(self):
        self.connection.close()
