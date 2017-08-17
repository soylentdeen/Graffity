import sys

sys.path.append('../')

import Graffity
from matplotlib import pyplot
from matplotlib.backends.backend_pdf import PdfPages
import CIAO_DatabaseTools
import numpy

CIAO_keywords = ['STREHL', 'SEEING', 'ALT', 'AZ']

CIAO_DB = CIAO_DatabaseTools.CIAO_Database()

CIAO_values = CIAO_DB.query(keywords=CIAO_keywords, timeOfDay='NIGHT',
        startTime='2017-08-08 00:00:00')

AVCModes = {}

frequencies = [24.5, 48.0, 74.0, 76.0, 78.0, 97.0, 99.0] 

Power={}
Alt = {}
Az = {}
Date = {}

#symbols = ['o', 'v', '^', '*', '<']
colors = [numpy.random.random(3,) for i in frequencies]
#colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'b', 'g', 'r']

#"""
for CIAO_ID in [1,2, 3, 4]:
    print CIAO_ID
    Power[CIAO_ID] = {}
    Alt[CIAO_ID] = {}
    Az[CIAO_ID] = {}
    Date[CIAO_ID] = {}
    for record in CIAO_values[CIAO_ID]:
        dl = Graffity.DataLogger(directory=record[-3])
        print record[-3]
        dl.loadData()
        dl.computeStrehl()
        dl.measureVibs(frequencies=frequencies, modes='AVC')
        for key in dl.vibPower.keys():
            if not(key in Power[CIAO_ID].keys()):
                Power[CIAO_ID][key] = {}
                Alt[CIAO_ID][key] = {}
                Az[CIAO_ID][key] = {}
                Date[CIAO_ID][key] = {}
            for f in dl.vibPower[key]["CommPower"].keys():
                if not(f in Power[CIAO_ID][key].keys()):
                    Power[CIAO_ID][key][f] = []
                    Alt[CIAO_ID][key][f] = []
                    Az[CIAO_ID][key][f] = []
                    Date[CIAO_ID][key][f] = []
                Power[CIAO_ID][key][f].append(dl.vibPower[key]["CommPower"][f])
                Alt[CIAO_ID][key][f].append(record[2])
                Az[CIAO_ID][key][f].append(record[3])
                Date[CIAO_ID][key][f].append(record[-4])

        del(dl)

Azpdf = PdfPages('Azimuth.pdf')
Altpdf = PdfPages('Altitude.pdf')
Timepdf = PdfPages('Time.pdf')

for mode, i in zip(Power[1].keys(), range(len(Power[1].keys()))):
    for pdf, ordinal in zip([Azpdf, Altpdf, Timepdf], [Az, Alt, Date]):
        fig = pyplot.figure(0, figsize=(5, 10))
        fig.clear()
        ax1 = fig.add_axes([0.1, 0.3, 0.4, 0.3])
        ax2 = fig.add_axes([0.5, 0.3, 0.4, 0.3])
        ax3 = fig.add_axes([0.1, 0.6, 0.4, 0.3])
        ax4 = fig.add_axes([0.5, 0.6, 0.4, 0.3])
        ax1.clear()
        ax2.clear()
        ax3.clear()
        ax4.clear()
        fig.suptitle('Mode '+str(i+1)+' Vibrations')
        ax1.text(0.05, 0.05, 'UT1', transform=ax1.transAxes)
        ax2.text(0.05, 0.05, 'UT2', transform=ax2.transAxes)
        ax3.text(0.05, 0.05, 'UT3', transform=ax3.transAxes)
        ax4.text(0.05, 0.05, 'UT4', transform=ax4.transAxes)
        for freq, j in zip(frequencies, range(len(frequencies))):
            for CIAO_ID, ax in zip([1, 2, 3, 4], [ax1, ax2, ax3, ax4]):
                ax.scatter(numpy.array(ordinal[CIAO_ID][mode][freq]),
                    numpy.log10(numpy.array(Power[CIAO_ID][mode][freq])),color=colors[j],
                    label="%.1f Hz" % freq, s=5)
        handles, labels = ax1.get_legend_handles_labels()
        fig.legend(handles, labels, scatterpoints=1, ncol=4, frameon=True, loc=3)
        pdf.savefig(fig)

Azpdf.close()
Altpdf.close()
Timepdf.close()
