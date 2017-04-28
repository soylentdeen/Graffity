import numpy
import scipy
from matplotlib import pyplot

def fitStrehl(ciao, iris):
    A = []
    for c in ciao:
        A.append([0.0, c])

    A = numpy.matrix(A)
    x = A.I.dot(iris)

    answer = A.dot(x.T).tolist()

    return numpy.array(answer)[:,0], numpy.array(x.tolist())

df = open('/home/cdeen/Data/CIAO/JanComm/IRIS_Data.dat', 'r')

Seeing = []
Strehl = []
Iris = []
Flux = []
Alt = []
Az = []
Ciao = []

c = 0
colors = numpy.array(['g', 'r', 'b', 'm'])
for line in df.readlines():
    if line[0] != 'C':
        if line[0] != '#':
            l = line.split('|')
            Seeing.append(numpy.float32(l[1].strip()))
            Strehl.append(numpy.float32(l[2].strip()))
            Iris.append(numpy.float32(l[6].strip()))
            Flux.append(numpy.float32(l[5].strip()))
            Alt.append(numpy.float32(l[3].strip()))
            Az.append(numpy.float32(l[4].strip()))
            Ciao.append(colors[c-1])
    else:
        c += 1


Seeing = numpy.array(Seeing)
Strehl = numpy.array(Strehl)
Iris = numpy.array(Iris)
Flux = numpy.array(Flux)
Alt = numpy.array(Alt)
Az = numpy.array(Az)
Ciao = numpy.array(Ciao)

fit1, f1 = fitStrehl(Strehl[Ciao=='g'], Iris[Ciao=='g'])
fit2, f2 = fitStrehl(Strehl[Ciao=='r'], Iris[Ciao=='r'])
fit3, f3 = fitStrehl(Strehl[Ciao=='b'], Iris[Ciao=='b'])
fit4, f4 = fitStrehl(Strehl[Ciao=='m'], Iris[Ciao=='m'])
fit5, f5 = fitStrehl(Strehl, Iris)

fig = pyplot.figure(0)
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
ax.clear()

ax.scatter(Strehl[Ciao == 'g'], Iris[Ciao== 'g'], color = Ciao[Ciao== 'g'], label='UT1')
ax.scatter(Strehl[Ciao == 'r'], Iris[Ciao== 'r'], color = Ciao[Ciao== 'r'], label='UT2')
ax.scatter(Strehl[Ciao == 'b'], Iris[Ciao== 'b'], color = Ciao[Ciao== 'b'], label='UT3')
ax.scatter(Strehl[Ciao == 'm'], Iris[Ciao== 'm'], color = Ciao[Ciao== 'm'], label='UT4')
#ax.scatter(Strehl[Ciao=='g'], fit1, c='g', marker = '+')
#ax.scatter(Strehl[Ciao=='r'], fit2, c='r', marker = '+')
#ax.scatter(Strehl[Ciao=='b'], fit3, c='b', marker = '+')
#ax.scatter(Strehl[Ciao=='m'], fit4, c='m', marker = '+')
#ax.scatter(Strehl, fit5, c='k', marker='^')
ax.plot([0.0, 1.0], [0.0, 0.716], label = 'Fit')
ax.legend(loc=2, scatterpoints = 1)
ax.set_xbound(0.0, 1.0)
ax.set_ybound(0.0, 0.8)
ax.set_xlabel("CIAO caculated Strehl ratio")
ax.set_ylabel("IRIS measured Strehl ratio")
ax.set_title("Calibration of CIAO Calculated Strehl Ratios")
ax.set_aspect('auto')
#ax.scatter(Seeing[Flux > 2000], Iris[Flux > 2000], color = Ciao[Flux > 2000])
#ax.scatter(Seeing[Flux > 1000], 0.66*Strehl[Flux > 1000], color = Ciao[Flux > 1000])
#ax.scatter(Seeing, Iris, color = Ciao)
#ax.scatter(Flux, Strehl, color = Ciao)
fig.show()
fig.savefig("StrehlCalibration.png")
