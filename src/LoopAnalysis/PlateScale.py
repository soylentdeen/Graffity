import Graffity
import numpy
import matplotlib.pyplot as pyplot
import scipy
import pyfits
import sys
import glob
from scipy import optimize

def fitCircle(x, y, theta):
    errfunc_X = lambda p, th : p[0] + p[2]*numpy.cos(th+p[4])
    errfunc_Y = lambda p, th : p[1] + p[3]*numpy.sin(th+p[4])
    #fitfunc = lambda p, x, y, th : numpy.abs(errfunc_X(p, th) - x) + numpy.abs(errfunc_Y(p, th) - y)
    fitfunc = lambda p, x, y, th : numpy.abs(errfunc_X(p, th) - x) + numpy.abs(errfunc_Y(p, th) - y)
    guessX = numpy.mean(x)
    guessY = numpy.mean(y)
    r = numpy.std(x)
    params = [guessX, guessY, r, r, 0.01]
    pfit = optimize.leastsq(fitfunc, params, args = (x, y, theta))
    return pfit[0]

datadir = '/home/deen/Data/GRAVITY/MATLAB_DATA/2016-03-12/DM_CONJ-200255/'

loops = glob.glob(datadir+'CIAO_LOOP*.fits')
loops.sort()

pyplot.rcParams['font.size'] = 18
fig = pyplot.figure(0)
fig.clear()
ax = fig.add_axes([0.15, 0.15, 0.8, 0.8])

tip = []
tilt = []
x = []
y = []

npts = len(loops)-1
delta = 2.0*numpy.pi/npts
theta = []

for loop in loops:
    if len(x) == 0:
        x.append(0.0)
        y.append(0.0)
    else:
        i = len(x)-1
        x.append(numpy.cos(i*delta))
        y.append(numpy.sin(i*delta))
        theta.append(i*delta)

    data = pyfits.getdata(loop)
    
    TTM = data.field("ITTM_Positions")

    tip.append(numpy.median(TTM[:,0]))
    tilt.append(numpy.median(TTM[:,1]))

circle = fitCircle(tip[1:], tilt[1:], theta)

centx = circle[0]
centy = circle[1]
rx = circle[2]
ry = circle[3]
points_x = []
points_y = []

newTheta = numpy.linspace(0, numpy.pi*2.0)

for th in newTheta:
    points_x.append(centx + rx*numpy.cos(th))
    points_y.append(centy + ry*numpy.sin(th))


ax.scatter(tip, tilt, s=30)
ax.plot(points_x, points_y, color = 'r', lw=3.0)

plateScale = numpy.mean([numpy.abs(rx), numpy.abs(ry)])

ax.text(0.2, 0.9, 'Plate Scale: %.3f TTMU/pix' % plateScale, 
        fontsize=20, transform=ax.transAxes)

ax.set_ylabel("Tilt Position")
ax.set_xlabel("Tip Position")

ax.set_aspect('equal')
#ax.scatter(x, y)
#ax.set_xbound(-0.1, 0.1)
#ax.set_ybound(-0.1, 0.1)
fig.show()
fig.savefig("PlateScale.png")
