import sys

sys.path.append('../')

import numpy
import Graffity
import CIAO_DatabaseTools
import astropy.time as aptime
from matplotlib import pyplot
import colorsys

def getFreqs():
    while True:
        retval = []
        enteredText = raw_input("Enter a comma separated list of frequencies: ")
        try:
            for val in enteredText.split(','):
                retval.append(float(val.strip()))
            break
        except:
            pass
    return retval


def getModes():
    while True:
        enteredText = raw_input("Which modes to investigate? AVC or ALL? : ")
        if enteredText == 'AVC':
            return 'AVC'
        if enteredText == 'ALL':
            return 'ALL'

def getDataLoggers(DB, GravityVals, startTime):
    order = numpy.argsort(GravityVals[:,-2])
    GravityVals = GravityVals[order]
    i = 1
    for record in GravityVals:
        print("%03d | %s" % (i,aptime.Time(float(record[-2]), format='mjd').iso))
        i += 1

    index = int(raw_input("Enter desired index :")) - 1
    FT_OPDS = Graffity.GRAVITY_Data(GravityVals[index][-1])
    FT_OPDS.DualSciP2VM.computeOPDPeriodograms()
    return FT_OPDS
    """
    print asdf
    freqs = getFreqs()
    Modes = getModes()

    CIAOVals = DB.query(keywords=['ALT', 'AZ', 'STREHL'], timeOfDay='NIGHT', startTime=startTime)

    DataLoggers = {}
    for UT in [1, 2, 3, 4]:
        closest = numpy.argsort(numpy.abs(CIAOVals[UT][:,-4]
         - float(GravityVals[index,-2])))[0]

        DataLoggers[UT] = Graffity.DataLogger(directory=CIAOVals[UT][closest,-3])
        DataLoggers[UT].loadData()
        DataLoggers[UT].computeStrehl()
        DataLoggers[UT].measureVibs(frequencies=freqs, modes=Modes)

    return DataLoggers
    """

fig = pyplot.figure(0, figsize=(8.0, 10.0), frameon=False)
fig.clear()
ax1 = fig.add_axes([0.1, 0.2, 0.4, 0.3])
ax2 = fig.add_axes([0.1, 0.5, 0.4, 0.4], sharex=ax1)
ax3 = fig.add_axes([0.5, 0.2, 0.4, 0.3], sharex=ax1)
ax3.yaxis.tick_right()
ax4 = fig.add_axes([0.5, 0.5, 0.4, 0.4], sharex=ax1)
ax4.yaxis.tick_right()

GDB = CIAO_DatabaseTools.GRAVITY_Database()
CDB = CIAO_DatabaseTools.CIAO_Database()

startTime = '2017-08-08 00:00:00'

GravityVals = GDB.query(keywords = [], timeOfDay='NIGHT', startTime=startTime)

CIAO = getDataLoggers(CDB, GravityVals, startTime)
print asdf
hsv = [(numpy.random.uniform(low=0.0, high=1),
           numpy.random.uniform(low=0.2, high=1),
           numpy.random.uniform(low=0.9, high=1)) for i in
           range(99)]
colors = []
for h in hsv:
    colors.append(colorsys.hsv_to_rgb(h[0], h[1], h[2]))

handles = []
labels = []
for CIAO_ID, ax in zip([1, 2, 3, 4], [ax1, ax2, ax3, ax4]):
    DL = CIAO[CIAO_ID]
    DL.pupilIllumination(ax)
    Scale = 1e9*numpy.sqrt(DL.ZPowerdFrequencies)
    """
    for mode in DL.vibPower.keys():
        f = []
        p = []
        for peak in DL.vibPower[mode]['CommPower'].iteritems():
            f.append(peak[0])
            #p.append(numpy.log10(peak[1]))
            #p.append(numpy.sqrt(peak[1])*Scale)
            p.append(peak[1])
        #ax.plot(DL.ZPowerFrequencies, numpy.log10(DL.ZPowerCommands[mode,:]), color =
        #        colors[mode])
        #ax.plot(DL.ZPowerFrequencies, numpy.sqrt(DL.ZPowerCommands[mode,:])*Scale, color =
        #        colors[mode])
        ax.plot(DL.ZPowerFrequencies, DL.ZPowerCommands[mode,:], color =
                colors[mode])
        ax.scatter(numpy.array(f), numpy.array(p), color=colors[mode],
                    label='Mode %d' % mode)
        #for freq in f:
        #    ax.plot([f, f], [0, 100], color='k', lw=0.1)
"""
#handles, labels = ax.get_legend_handles_labels()
#ax1.set_ybound(0, 20)
#ax2.set_ybound(0, 20)
#ax3.set_ybound(0, 20)
#ax4.set_ybound(0, 20)
ax1.set_xbound(0, 160)
ax2.set_xbound(0, 160)
ax3.set_xbound(0, 160)
ax4.set_xbound(0, 160)
ax1.text(0.05, 0.05, 'UT1', transform=ax1.transAxes)
ax2.text(0.05, 0.05, 'UT2', transform=ax2.transAxes)
ax3.text(0.05, 0.05, 'UT3', transform=ax3.transAxes)
ax4.text(0.05, 0.05, 'UT4', transform=ax4.transAxes)
#ax2.xaxis.set_ticklabels([])
#ax4.xaxis.set_ticklabels([])

#fig.legend(handles, labels, ncol=4, loc=3, scatterpoints=1)
#"""
ax1.set_xlabel("Frequency (Hz)")
ax1.set_ylabel("FFT of Edge Subaperture Illumination")
fig.show()  
fig.savefig("PupilIllumination.png")
