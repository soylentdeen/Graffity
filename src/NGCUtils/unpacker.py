import pyfits
import sys
import numpy
import matplotlib.pyplot as pyplot

fig = pyplot.figure(0)
fig.clear()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

df = sys.argv[1]

cube = pyfits.open(df)

header = cube[0].header

medians = []
means = []

for frame in cube[0].data:
    medians.append(numpy.median(frame))
    means.append(numpy.mean(frame))


ax.plot(medians)
ax.plot(means)
fig.show()
