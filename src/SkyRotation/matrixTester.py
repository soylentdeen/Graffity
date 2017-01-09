import numpy
import pyfits
import matplotlib.pyplot as pyplot

stat = pyfits.getdata('/home/deen/Data/GRAVITY/SPARTA_Data/Matrices/UT1/IMs/STAT/20160315T193144/RTC.STAT.HOIM.fits', ignore_missing_end=True)
sky = pyfits.getdata('/home/deen/Data/GRAVITY/SPARTA_Data/Matrices/UT1/IMs/SKY/20160315T193144/RTC.SKY.HOIM.fits', ignore_missing_end=True)

stat_0 = stat[:,0]
diff = []
for sky_0 in sky.T:
    diff.append(numpy.sum(numpy.abs(stat_0-sky_0)))

fig = pyplot.figure(0)
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

ax.plot(diff)

fig.show()



