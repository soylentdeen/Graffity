import scipy
import scipy.linalg as linalg
import numpy
import matplotlib.pyplot as pyplot

datafiles = ['UT3.dat', 'UT1.dat']

fig = pyplot.figure(0)
fig.clear()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

for datafile in datafiles:
    VCM = []
    FLDL = []
    data = open(datafile, 'r')
    for line in data.readlines():
        if line[0] != '#':
             vals = numpy.array([x.strip() for x in line.split(',')], dtype=numpy.float32)
             VCM.append(numpy.array([vals[0], vals[1]]))
             FLDL.append(numpy.array([vals[2], vals[3]]))

    VCM = numpy.array(VCM)
    FLDL = numpy.array(FLDL)

    #meanVCM = numpy.mean(VCM, axis=0)
    #meanFLDL = numpy.mean(FLDL, axis=0)
    VCM = VCM# - meanVCM
    FLDL = FLDL# - meanFLDL

    #x = linalg.pinv(VCM).dot(FLDL)
    x = scipy.matrix([[1.0, x[0], x[1]] for x in VCM])

    mat = linalg.pinv(x.T.dot(x)).dot(x.T).dot(FLDL)
    answer = mat.T.dot(x.T).T
    print mat.T
    #print x, datafile

    #print (VCM).dot(x) - FLDL

    for v, f, a in zip(VCM, FLDL, answer[:,:]):
        ax.plot([v[0], f[0]], [v[1], f[1]], color = 'b')
        ax.plot([v[0], a.flat[0]], [v[1], a.flat[1]], color= 'r')

    fig.show()
    raw_input()
    ax.clear()

