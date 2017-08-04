import scipy
import numpy
import sys
from astropy import time as aptime

sys.path.append('../')

import CIAO_DatabaseTools
import Graffity
import tkinter as tk
import ttk
from matplotlib import pyplot
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import matplotlib.backend_bases
from matplotlib.figure import Figure


class Browser( object ):
    def __init__(self, GRAVITY_values, CIAO_keywords, GRAVITY_keywords, CIAO_DB, startTime):
         
        self.root = tk.Tk()
        self.root.title("CIAO/GRAVITY Observation Browser")

        self.GRAVITY_values = GRAVITY_values
        self.CIAO_keywords = CIAO_keywords
        self.GRAVITY_keywords = GRAVITY_keywords
        self.CIAO_DB = CIAO_DB
        self.startTime = startTime


        self.Obs_Canvas = tk.Canvas(self.root, borderwidth=0, background='#ffffff',
                width=300, height=900)
        self.Plot_Canvas = tk.Canvas(self.root, borderwidth=0, background='#ffffff',
                width=800, height=900)
        self.frame1 = tk.Frame(self.Obs_Canvas, background='#ffffff')
        self.frame2 = tk.Frame(self.Plot_Canvas, background='#ffffff')
        self.frame3 = tk.Frame(self.Plot_Canvas, background='#ffffff')
        self.vsb = tk.Scrollbar(self.root, orient='vertical', command=self.Obs_Canvas.yview)
        self.Obs_Canvas.configure(yscrollcommand=self.vsb.set)

        self.Obs_Canvas.pack(side='left', fill='both', expand=False)
        self.vsb.pack(side='left', fill='y')
        self.Obs_Canvas.create_window((4,4), window=self.frame1, anchor='nw',
                tags='self.frame1')
        self.Plot_Canvas.pack(side='left', expand=True)
        self.Plot_Canvas.create_window((4, 4), window=self.frame2, anchor='nw',
                tags='self.frame2')
        self.Plot_Canvas.create_window((4, 400), window=self.frame3, anchor='nw',
                tags='self.frame3')
        self.colors = {1:'g', 2:'b', 3:'r', 4:'k'}

        self.frame1.bind("<Configure>", self.onFrameConfigure)

        self.populate()

        self.root.mainloop()
    
    def populate(self):
        tk.Label(self.frame1, text="GRAVITY Observations", width = 30, borderwidth="1",
                relief="solid").grid(row=0, column = 0)
        i = 0
        order = numpy.argsort(GRAVITY_values[:, -2])
        self.ObsCheckboxVars = []
        self.ObsTime = []
        for val in self.GRAVITY_values[order]:
            self.ObsCheckboxVars.append(tk.IntVar())
            self.ObsTime.append(float(val[-2]))
            tk.Checkbutton(self.frame1, text="%03d | %s " % (i,
                aptime.Time(float(val[-2]), format='mjd').iso),
                variable=self.ObsCheckboxVars[-1]).grid(row=i+1, column=0)
            i += 1

        tk.Label(self.frame2, text="CIAO Observations", width=30, borderwidth="1", 
                relief="solid").grid(row=0, column=0, columnspan=6)
        self.UTCheckboxVars = {}
        for i, c in zip([1, 2, 3, 4], ['#0a9d1f', '#0e12bc', '#ff0000', '#000000']):
            self.UTCheckboxVars[i] = tk.IntVar()
            cb = tk.Checkbutton(self.frame2, text="UT%d" % i, bg=c,
                    variable=self.UTCheckboxVars[i])
            cb.select()
            cb.grid(row=1, column=i-1)
        tk.Button(self.frame2, text="Quit", command=self.root.quit).grid(row=1,
                column=4)
        tk.Button(self.frame2, text="Show", command=self.showData).grid(row=1,
                column = 5)

        self.CIAO_nb = ttk.Notebook(self.frame2)
        self.CIAO_axes = {}
        self.CIAO_tabs = {}
        for kw in self.CIAO_keywords:
            print kw
            f = Figure(figsize=(5,4))
            self.CIAO_axes[kw] = f.add_axes([0.1, 0.1, 0.8, 0.8])
            self.CIAO_axes[kw].set_title(kw)
            self.CIAO_tabs[kw] = ttk.Frame(self.CIAO_nb, width=600, height=600)
            canvas = FigureCanvasTkAgg(f, master=self.CIAO_tabs[kw])
            canvas.show()
            canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
            self.CIAO_nb.add(self.CIAO_tabs[kw], text=kw)

        self.CIAO_nb.grid(row=2,column=0, columnspan=6)

        tk.Label(self.frame3, text="GRAVITY Observation values", width=30,
                borderwidth="1", relief="solid").grid(row=0, column=0)

        self.GRAVITY_nb = ttk.Notebook(self.frame3)
        self.GRAVITY_axes = {}
        self.GRAVITY_tabs = {}

        for kw in self.GRAVITY_keywords:
            print kw
            f = Figure(figsize=(5, 4))
            self.GRAVITY_axes[kw] = f.add_axes([0.1, 0.1, 0.8, 0.8])
            self.GRAVITY_axes[kw].set_title(kw)
            self.GRAVITY_tabs[kw] = ttk.Frame(self.GRAVITY_nb, width=600,
                    height=600)
            canvas = FigureCanvasTkAgg(f, master=self.GRAVITY_tabs[kw])
            canvas.show()
            canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
            self.GRAVITY_nb.add(self.GRAVITY_tabs[kw], text=kw)

        self.GRAVITY_nb.grid(row=1,column=0)
        
    def onFrameConfigure(self, event):
        self.Obs_Canvas.configure(scrollregion=self.Obs_Canvas.bbox("all"))


    def showData(self):
        times = {}
        i = 0
        for time, cb in zip(self.ObsTime, self.ObsCheckboxVars):
            if cb.get() == 1:
                times[i] = time
                i += 1

        CIAO_values = self.CIAO_DB.query(keywords=self.CIAO_keywords, timeOfDay='NIGHT',
                startTime=self.startTime)

        DataLoggers = {1:[], 2:[], 3:[], 4:[]}
        CIAO_Values = {}
        for key in self.CIAO_keywords:
            CIAO_Values[key] = {1:[], 2:[], 3:[], 4:[]}
        for i in times.keys():
            for UT in self.UTCheckboxVars.keys():
                if self.UTCheckboxVars[UT].get() == 1:
                    timeStamp = numpy.argsort(numpy.abs(numpy.array(CIAO_values[UT][:,-4],
                        dtype=numpy.float32) - float(times[i])))[0]
                    timeDistance = (float(CIAO_values[UT][timeStamp, -4]) - float(times[i]))*24*3600
                    if timeDistance < 90:
                        print("Distance for UT %d : %.3f seconds" % (UT, timeDistance))
                        DataLoggers[UT].append(Graffity.DataLogger(directory=CIAO_values[UT][timeStamp,-3]))
                        DataLoggers[UT][-1].loadData()
                        for key, j in zip(CIAO_keywords, range(len(CIAO_keywords))):
                            try:
                                CIAO_Values[key][UT].append(float(CIAO_values[UT][timeStamp,j]))
                            except:
                                CIAO_Values[key][UT].append(0.0)
                    else:
                        print("Error!  Datalogger within 30 seconds does not exist for this observation!")

        for kw in self.CIAO_keywords:
            self.CIAO_axes[kw].clear()
            for UT in CIAO_Values[kw].keys():
                if self.UTCheckboxVars[UT].get() == 1:
                    self.CIAO_axes[kw].scatter(times.keys(), CIAO_Values[kw][UT], color=self.colors[UT])
            self.CIAO_axes[kw].figure.canvas.draw()

        for kw in self.GRAVITY_keywords:
            self.GRAVITY_axes[kw].clear()
            for UT in GRAVITY_Values[kw].keys():
                if self.UTCheckboxVars[UT].get() == 1:
                    self.GRAVITY_axes[kw].scatter(times.keys(), GRAVITY_Values[kw][UT], color=self.colors[UT])
            self.CIAO_axes[kw].figure.canvas.draw()



if __name__ == "__main__":
    CIAO_keywords= ['STREHL', 'SEEING', 'ASM_SEEING', 'M10_POSANG',
                  'WINDDIR', 'WINDSP', 'PRLTIC', 'TIP_RESIDUALS',
                  'TILT_RESIDUALS', 'ALT', 'AZ']
    GRAVITY_keywords=['ACQCAM_1_STREHL', 'ACQCAM_2_STREHL', 'ACQCAM_3_STREHL',
    'ACQCAM_4_STREHL']
    CIAO_DB = CIAO_DatabaseTools.CIAO_Database()
    GRAVITY_DB = CIAO_DatabaseTools.GRAVITY_Database()
    
    startTime = '2017-08-01 00:00:00'
    
    GRAVITY_values = GRAVITY_DB.query(keywords = [], timeOfDay='NIGHT',
            startTime=startTime)
    Browser(GRAVITY_values, CIAO_keywords, GRAVITY_keywords, CIAO_DB, startTime)

