import scipy
import numpy
import sys
from astropy import time as aptime

sys.path.append('../')

import CIAO_DatabaseTools
import Graffity
import tkinter

def getGRAVITY_OBS(GRAVITY_values, frame):
    tkinter.Label(frame, text="GRAVITY Observations", width = 3, borderwidth="1",
            relief="solid").grid(row=0, column = 0)
    """"
    print("GRAVITY Observations")
    i = 0
    print("i | Filename | Strehl | Seeing")
    order = numpy.argsort(GRAVITY_values[:, -2])
    for val in GRAVITY_values[order]:
        print("%03d | %s " % (i, aptime.Time(float(val[-2]), format='mjd').iso))
        i += 1
        

    choice = raw_input("Enter comma separated choices ('END' to quit): ")
    try:
        choices = [int(r.strip()) for r in choice.split(',')]
        retval = []
        for c in choices:
            retval.append(GRAVITY_values[order][c,-2])
        return retval
    except:
        return "END"
    """

def getDataLoggers(GRAVITY_OBS, CIAO_DB, UTS=[1,2,3,4]):
    keywords= ['STREHL', 'SEEING', 'ASM_SEEING', 'M10_POSANG',
                  'WINDDIR', 'WINDSP', 'PRLTIC', 'TIP_RESIDUALS',
                  'TILT_RESIDUALS', 'ALT', 'AZ']
    CIAO_values = CIAO_DB.query(keywords= keywords, timeOfDay='NIGHT',
            startTime=startTime)
    DataLoggers = {1:[], 2:[], 3:[], 4:[]}
    Values = {}
    for key in keywords:
        Values[key] = {1:[], 2:[], 3:[], 4:[]}
    for GRAVOBS in GRAVITY_OBS:
        for UT in UTS:
            timeStamp = numpy.argsort(numpy.abs(numpy.array(CIAO_values[UT][:,-4],
                dtype=numpy.float32) - float(GRAVOBS)))[0]
            timeDistance = (float(CIAO_values[UT][timeStamp, -4]) - float(GRAVOBS))*24*3600
            if timeDistance < 90:
                print("Distance for UT %d : %.3f seconds" % (UT, timeDistance))
                DataLoggers[UT].append(Graffity.DataLogger(directory=CIAO_values[UT][timeStamp,-3]))
                DataLoggers[UT][-1].loadData()
                for key, i in zip(keywords, range(len(keywords))):
                    try:
                        Values[key][UT].append(float(CIAO_values[UT][timeStamp,i]))
                    except:
                        Values[key][UT].append(0.0)
            else:
                print("Error!  Datalogger within 30 seconds does not exist for this observation!")
    return Values, DataLoggers

def onFrameConfigure(canvas):
    canvas.configure(scrollregion=canvas.bbox("all"))

Root = tkinter.Tk()

CIAO_DB = CIAO_DatabaseTools.CIAO_Database()
GRAVITY_DB = CIAO_DatabaseTools.GRAVITY_Database()

UTS = [2]

startTime = '2017-07-01 00:00:00'

GRAVITY_values = GRAVITY_DB.query(keywords = [], timeOfDay='NIGHT',
        startTime=startTime)

canvas = tkinter.Canvas(Root, borderwidth=0, background='#ffffff')
frame = tkinter.Frame(canvas, background='#ffffff')
vsb = tkinter.Scrollbar(Root, orient='vertical', command=canvas.yview)
canvas.configure(yscrollcommand=vsb.set)

vsb.pack(side='right', fill='y')
canvas.pack(side='left', fill='both', expand=True)
canvas.create_window((4,4), window=frame, anchor='nw')

frame.bind("<Configure>", lambda event, canvas=canvas: onFrameConfigure(canvas))

#getGRAVITY_OBS(GRAVITY_values, frame)

Root.mainloop(frame)

"""
while True:
    GRAVITY_OBS = getGRAVITY_OBS(GRAVITY_values)
    if GRAVITY_OBS == "END":
        break
    CIAO_OBS, CIAO_DLS = getDataLoggers(GRAVITY_OBS, CIAO_DB)

    print asdf


"""
