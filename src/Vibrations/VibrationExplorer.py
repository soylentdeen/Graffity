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

def getDataLoggers(DB, GravityVals, startTime, ax=None):
    order = numpy.argsort(GravityVals[:,-2])
    GravityVals = GravityVals[order]
    i = 1
    for record in GravityVals:
        print("%03d | %s" % (i,aptime.Time(float(record[-2]), format='mjd').iso))
        i += 1

    index = int(raw_input("Enter desired index :")) - 1
    FTData = Graffity.GRAVITY_Data(GravityVals[index][-1])
    FTData.DualSciP2VM.computeOPDPeriodograms()
    VibrationPeaks = FTData.DualSciP2VM.findVibrationPeaks()

    FTData.computeACQCAMStrehl()

    #freqs = getFreqs()
    #Modes = getModes()

    CIAOVals = DB.query(keywords=['ALT', 'AZ', 'STREHL'], timeOfDay='NIGHT', startTime=startTime)

    DataLoggers = {}
    for UT in [1, 2, 3, 4]:
        closest = numpy.argsort(numpy.abs(CIAOVals[UT][:,-4]
         - float(GravityVals[index,-2])))[0]

        DataLoggers[UT] = Graffity.DataLogger(directory=CIAOVals[UT][closest,-3])
        DataLoggers[UT].loadData()
        DataLoggers[UT].computeStrehl()
        freqs = extractBCIFreqs(VibrationPeaks, UT)
        DataLoggers[UT].measureVibs(frequencies=freqs, modes='AVC')

    return DataLoggers, VibrationPeaks

def extractBCIFreqs(VibrationPeaks, UT):
    freqs = []
    baselines = {0:[4,3], 1:[4, 2], 2:[4, 1], 3:[3, 2], 4:[3, 1], 5:[2, 1]}
    for bl in baselines.keys():
        if UT in baselines[bl]:
            for f in VibrationPeaks[bl]['freqs']:
                freqs.append(f)
    return numpy.array(freqs)


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

startTime = '2017-08-10 00:00:00'

GravityVals = GDB.query(keywords = [], timeOfDay='NIGHT', startTime=startTime)

#ax1.set_xscale('log')
#ax1.set_yscale('log')
CIAO, Vibrations = getDataLoggers(CDB, GravityVals, startTime, ax=ax1)

hsv = [(numpy.random.uniform(low=0.0, high=1),
           numpy.random.uniform(low=0.2, high=1),
           numpy.random.uniform(low=0.9, high=1)) for i in
           range(99)]
colors = []
for h in hsv:
    colors.append(colorsys.hsv_to_rgb(h[0], h[1], h[2]))

handles = numpy.array([])
labels = numpy.array([])
baselines = {0:[4,3], 1:[4, 2], 2:[4, 1], 3:[3, 2], 4:[3, 1], 5:[2, 1]}
colors = {0:'y', 1:'g', 2:'r', 3:'c', 4:'m', 5:'k'}
for CIAO_ID, ax in zip([1, 2, 3, 4], [ax1, ax2, ax3, ax4]):
    DL = CIAO[CIAO_ID]
    for mode in DL.vibPower.keys():
        BCIVibs = {}
        for bl in baselines.keys():
            if CIAO_ID in baselines[bl]:
                label = "UT%dUT%d" % (baselines[bl][0], baselines[bl][1])
                BCIVibs[label] = {'index':bl, 'power':[]}
        f = []
        p = []
        for peak in DL.vibPower[mode]['CommPower'].iteritems():
            if peak[1] > 0:
                f.append(peak[0])
                p.append(numpy.log10(peak[1]))
                for label in BCIVibs.keys():
                    if not( f[-1] in Vibrations[BCIVibs[label]['index']]['freqs']):
                        BCIVibs[label]['power'].append(0.0)
                    else:
                        for i, freq in enumerate(Vibrations[BCIVibs[label]['index']]['freqs']):
                            if freq == f[-1]:
                                BCIVibs[label]['power'].append(Vibrations[BCIVibs[label]['index']]['power'][i])
        
        #ax.plot(DL.ZPowerFrequencies, numpy.log10(DL.ZPowerCommands[mode,:]), color =
        #        colors[mode])
        f = numpy.array(f)
        p = numpy.array(p)
        ax.scatter(numpy.log10(f), p, color='b')
        for bl in BCIVibs.keys():
            BCIVibs[bl]['power'] = numpy.array(BCIVibs[bl]['power'])
            nonzero = BCIVibs[bl]['power'] > 0.0
            ax.scatter(numpy.log10(f[nonzero]), numpy.log10(BCIVibs[bl]['power'][nonzero]),
                    label=bl, color = colors[BCIVibs[bl]['index']])
        #ax.scatter(numpy.array(f), numpy.array(p), color=colors[mode],
        #            label='Mode %d' % mode)
        
        h, l = ax.get_legend_handles_labels()
        handles=numpy.append(handles, numpy.array(h))
        labels =numpy.append(labels, numpy.array(l))


#ax1.set_ybound(0, 20)
#ax2.set_ybound(0, 20)
#ax3.set_ybound(0, 20)
#ax4.set_ybound(0, 20)
#ax1.set_xbound(0, 160)
#ax2.set_xbound(0, 160)
#ax3.set_xbound(0, 160)
#ax4.set_xbound(0, 160)
#ax2.xaxis.set_ticklabels([])
#ax4.xaxis.set_ticklabels([])

junk, indices = numpy.unique(labels, return_index=True)

fig.legend(handles[indices], labels[indices], ncol=4, loc=3, scatterpoints=1)
fig.show()  
#"""
