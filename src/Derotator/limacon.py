import scipy
import numpy
import matplotlib.pyplot as pyplot

fig1 = pyplot.figure(1)
ax1 = fig1.add_axes([0.1, 0.1, 0.8, 0.8])
#ax1.clear()

theta = numpy.linspace(0, 2.0*numpy.pi)

p = numpy.array([0.0, 0.5, 2.0, 1.0, 0.0, 2.0])

C = 1.0 + 2.5j

A = -1.0 + 3.0j
B = 1.01 + 0.5j

ext = 0.00
internal = 0.000

I = C + A*numpy.exp((0.0+1.0j)*(2.0*theta + ext)) + B*numpy.exp((0.0+1.0j)*(theta+internal))

#x = p[0] + p[1] * numpy.cos(theta) + p[2]*numpy.cos(2.0*theta) - \
#        p[5]*numpy.sin(2.0*theta) - p[4]*numpy.sin(theta)
#y = p[3] + p[4] * numpy.sin(theta) + p[5]*numpy.sin(2.0*theta) - \

ax1.scatter(I.real, I.imag)

fig1.show()
