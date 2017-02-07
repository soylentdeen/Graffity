import scipy
import numpy
from scipy import linalg
import numpy
from matplotlib import pyplot
import glob
import os

def parseScriptName(name):
    extracted = name.split('.script')[0].split('_')

    value = float(extracted[1])
    if extracted[0] == 'V':
        return numpy.array([value, 0.0])/numpy.abs(value)
    else:
        return numpy.array([0.0, value])/numpy.abs(value)

def parseScriptFile(log):
    f = open(log, 'r')
    U = 0.0
    W = 0.0
    for line in f.readlines():
        split = line.split('mirrorId')
        if len(split) > 1:
            junk = split[1].split()
            U += float(junk[2])
            W += float(junk[4].strip("\'"))

    return numpy.array([U, W])/((U**2.0 + W**2.0)**0.5)

            

fig = pyplot.figure(0)
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

datadir = '/home/cdeen/Data/CIAO/UT1/VCM_Debugging/'

angle_dirs = [o for o in os.listdir(datadir) if os.path.isdir(os.path.join(datadir, o))]

angles = []
mats = []
X = []
b = []

for angle in angle_dirs:
    angles.append(float(angle))
    logs = glob.glob(os.path.join(datadir, angle)+'/*.script')
    for log in logs:
        X.append(numpy.append(parseScriptName(os.path.basename(log)), numpy.array([angles[-1]])))
        b.append(parseScriptFile(log))
X = numpy.array(X)
b = numpy.array(b)

#mat = linalg.pinv(b.T.dot(b)).dot(b.T).dot(X)
mat = linalg.pinv(X.T.dot(X)).dot(X.T).dot(b)


"""    
angles = numpy.array(angles)
order = numpy.argsort(angles)
angles = angles[order]
mats = numpy.array(mats)[order]

ax.clear()
ax.plot(angles, mats[:,0,0])
ax.plot(angles, mats[:,1,0])
ax.plot(angles, mats[:,0,1])
ax.plot(angles, mats[:,1,1])
fig.show()
"""
