import Graffity
import CIAO_DatabaseTools
from matplotlib import pyplot
import scipy
import numpy

db = CIAO_DatabaseTools.CIAO_Database()

UTS = [1, 2, 3, 4]

dataloggers = db.query(keywords = ["Strehl"], AVC_State='ON', timeOfDay="NIGHT", UTS=UTS, startTime='2017-03-25 00:00:00')
NCPA = {}
colors = ['b', 'g', 'r', 'c']

for UT in UTS:
    print UT
    NCPA[UT] = []
    for dl in dataloggers[UT]:
        DL = Graffity.DataLogger(directory = dl[-3])
        DL.loadData()
        NCPA[UT].append(DL.getRefSlopeZernikes())
        del(DL)
    NCPA[UT] = numpy.array(NCPA[UT])

fig = pyplot.figure(0)
fig.clear()
ax1 = fig.add_subplot(331)
ax1.clear()
ax2 = fig.add_subplot(332)
ax2.clear()
ax3 = fig.add_subplot(333)
ax3.clear()
ax4 = fig.add_subplot(334)
ax4.clear()
ax5 = fig.add_subplot(335)
ax5.clear()
ax6 = fig.add_subplot(336)
ax6.clear()
ax7 = fig.add_subplot(337)
ax7.clear()
ax8 = fig.add_subplot(338)
ax8.clear()
ax9 = fig.add_subplot(339)
ax9.clear()

for UT, color in zip(UTS, colors):
    order = numpy.argsort(dataloggers[UT][:,-4])
    for i, ax in zip(range(9), [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9]):
        ax.plot(dataloggers[UT][order,-4], NCPA[UT][order, i], color=color)
        ax.set_ybound(numpy.min(NCPA[UT][order, i]), numpy.max(NCPA[UT][order,i]))

fig.show()
