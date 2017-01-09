import scipy
import numpy
import matplotlib.pyplot as pyplot
import pyfits
from scipy import optimize

class postage_stamp( object ):
    def __init__(self, nx=8, ny=8):
        self.image = numpy.zeros((nx, ny), dtype=numpy.float32)
        self.xcenters = []
        self.ycenters = []
        self.fwhms = []

    def add(self, data):
        #fit = self.fit_gaussian(data)
        scaling_factor = numpy.max(data)
        self.image += data/scaling_factor

    def show(self):
        fig = pyplot.figure(0)
        fig.clear()
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        ax.imshow(self.image)
        fig.show()

def fit_line(x, y):
    fitfunc = lambda p, x : p[0]+(x*p[1])
    errfunc = lambda p, x, y: numpy.abs(fitfunc(p,x) - y)
    coeffs = [numpy.mean(y), 0.1]
    pfit = optimize.leastsq(errfunc, coeffs, args=(x,y) )

    #return numpy.arctan(numpy.abs(pfit[0][1]))*180.0/3.14159262
    return pfit

def sufficient_flux(x, y, img):
    median = 16*numpy.median(img)
    flux = numpy.sum(img[round(y-2):round(y+2),round(x-2):round(x+2)])
    #print x, y
    #print median, flux
    #raw_input()
    return flux > median

def calc_cog(img):
    nx = len(img)
    ny = len(img[0])
    valx = 0.0
    valy = 0.0
    tot = 0.0
    for i in range(nx):
        for j in range(ny):
            valx += j*img[j][i]
            valy += i*img[j][i]
            tot += img[j][i]


    return valx/tot, valy/tot

def extract(x, y, img):
    background = numpy.median(img)
    xcorner = round(x-3)
    ycorner = round(y-3)
    stamp = img[ycorner:ycorner+8,xcorner:xcorner+8].copy()
    stamp -= background

    order = numpy.argsort(numpy.ravel(stamp))
    stamp.ravel()[order[:-6]] = 0.0
    #xcog, ycog = calc_cog(stamp)
    xcoll = numpy.sum(stamp, axis=0)
    ycoll = numpy.sum(stamp, axis=1)
    xcog = numpy.sum(xcoll*numpy.arange(1,9))/numpy.sum(xcoll)-1.0
    ycog = numpy.sum(ycoll*numpy.arange(1,9))/numpy.sum(ycoll)-1.0
    """
    f1 = pyplot.figure(1)
    f1.clear()
    a1 = f1.add_axes([0.3, 0.1, 0.6, 0.6])
    a2 = f1.add_axes([0.1, 0.1, 0.2, 0.6])
    a3 = f1.add_axes([0.3, 0.7, 0.6, 0.2])
    a1.clear()
    a1.matshow(stamp)
    a1.scatter([4.0], [4.0], marker='o', s=10, color='r')
    a1.scatter([xcog], [ycog], marker='o')
    a2.plot(ycoll, range(len(stamp[0])))
    a2.plot(a2.get_xlim(), [ycog, ycog])
    a3.plot(range(len(stamp)), xcoll)
    a3.plot([xcog, xcog], a3.get_ylim())
    f1.show()
    print xcorner+xcog, ycorner+ycog
    raw_input()
    #"""
    return stamp, xcorner+xcog, ycorner+ycog

def subtract_gradient(img):
    background = img[:,0:70]
    gradient = background.sum(axis=1)/len(background[0,:])
    cleaned = img.copy()
    for i in range(len(img[0])):
        cleaned[:,i] -= gradient

    return cleaned



subaperture_positions_file = open('subaperture_positions_e.dat', 'r')
xpositions = numpy.array(subaperture_positions_file.readline().split(), 
        dtype=numpy.float32)
ypositions = numpy.array(subaperture_positions_file.readline().split(),
        dtype=numpy.float32)
subaperture_positions_file.close()

#datadir = '/home/deen/Data/GRAVITY/DetectorTests/MOVPE1/'
#df = datadir+'Light_MOVPE_DIT0.05_BIAS4.8_14-22.fits'
#datadir = '/home/deen/Data/GRAVITY/DetectorTests/MOVPE1/MOVPE_LLA#1_19.06.2015/'
#datadir = '/home/deen/Data/GRAVITY/DetectorTests/MOVPE1/Iteration2/'
datadir = '/home/deen/Data/GRAVITY/DetectorTests/MOVPE1/Iteration4/'
#df = datadir+'Light_MOVPE_DIT0.05_BIAS4.8_14-32.fits'
#df = datadir+'Light_MOVPE_DIT0.05_BIAS4.5_19-19.fits'
light_on_df = datadir+'0001_Light-on_Double_DIT0.05_Bias4.5.fits'
light_off_df = datadir+'0001_Light-off_Double_DIT0.05_Bias4.5.fits'
#df = 'clocking_image.fits'
data = pyfits.getdata(light_on_df)-pyfits.getdata(light_off_df)

#cleaned = subtract_gradient(data)
cleaned = data

stacked = postage_stamp()
xc = []
yc = []

for x in xpositions:
    for y in ypositions:
        if sufficient_flux(x, y, cleaned):
            extracted, cogx, cogy = extract(x, y, cleaned)
            stacked.add(extracted)
            xc.append(cogx)
            yc.append(cogy)
            #xc.append(x)
            #yc.append(y)

xc = numpy.array(xc)
yc = numpy.array(yc)
xangle = []
yangle = []

fig = pyplot.figure(0)
fig.clear()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

#ax.plot(cleaned.sum(axis=0))
#fig.show()
#raw_input()

meanx = numpy.tan(3.14159/180.0*0.1706)
meany = numpy.tan(3.15159/180.0*0.1633)

differences = []
for x in xpositions:
    col = scipy.where(abs(xc - x) < 3.0)
    if len(col[0]) > 5:
        fit = fit_line(yc[col], xc[col])
        ax.plot(fit[0][0]+fit[0][1]*yc[col], yc[col], color = 'r')
        ax.plot(fit[0][0]+meanx*yc[col], yc[col], color = 'k')
        ax.plot(numpy.ones(len(xc[col]))*numpy.mean(xc[col]), yc[col], color = 'b')
        #ax.scatter(xc[col], yc[col])
        differences.append(fit[0][0]+fit[0][1]*yc[col] - xc[col])
        #ax.plot(yc[col], fit[0][0]+fit[0][1]*yc[col] - xc[col])
        #numpy.arctan(numpy.abs(pfit[0][1]))*180.0/3.14159262
        xangle.append(numpy.arctan(numpy.abs(fit[0][1]))*180.0/3.14159)
        #print asdf

for y in ypositions:
    row = scipy.where(abs(yc - y) < 3.0)
    if len(row[0]) > 5:
        fit = fit_line(xc[row], yc[row])
        ax.plot(xc[row], fit[0][0]+fit[0][1]*xc[row], color = 'r')
        ax.plot(xc[row], fit[0][0]+meany*xc[row], color = 'k')
        ax.plot(xc[row], numpy.ones(len(yc[row]))*numpy.mean(yc[row]), color = 'b')
        #ax.scatter(xc[row], yc[row])
        #ax.plot(xc[row], fit[0][0]+fit[0][1]*xc[row] - yc[row])
        yangle.append(numpy.arctan(numpy.abs(fit[0][1]))*180.0/3.14159)
        #yangle.append(fit_line(xc[row], yc[row]))

print numpy.mean(xangle), numpy.std(xangle)
print numpy.mean(yangle), numpy.std(yangle)
#ax.plot(cleaned.sum(axis=0)/len(cleaned[:,0]))
#ax.clear()
#ax.matshow(cleaned)
ax.scatter(xc, yc, color = 'y', s=180)
fig.show()

fig.savefig("MOVPE_1_iteration4.png")
#stacked.show()

f2 = pyplot.figure(1)
ax2 = f2.add_axes([0.1, 0.1, 0.8, 0.8])
residual_x = xc % 1
residual_y = yc % 1

residual_x[residual_x > 0.5] -= 1.0
residual_y[residual_y > 0.5] -= 1.0

ax2.scatter(residual_x, residual_y)
print numpy.mean(residual_x)
print numpy.mean(residual_y)
f2.show()
f2.savefig("remainders_Iteration4.png")
