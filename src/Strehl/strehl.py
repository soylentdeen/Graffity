import scipy
import numpy
from matplotlib import pyplot

def generateOTF(sizeInPix = 20, lam = 2.0, dlam = 0.03, pscale = 17.0, M1=8.0, M2=1.3, nLambdaSteps = 9):
    lam = lam*1.0e-6
    dlam = dlam*1.0e-6
    pscale = pscale * 2.0 * numpy.pi/360.0/60.0**2.0/1.0e3
    fmax = M1 * pscale * sizeInPix / lam

    total = numpy.zeros([sizeInPix, sizeInPix])
    for k in numpy.arange(nLambdaSteps):
        l = lam - dlam*(k-nLambdaSteps/2.0)/nLambdaSteps
        fc = fmax*lam/l

        for i in range(sizeInPix):
            for j in range(sizeInPix):
                y = j - sizeInPix/2.0 + 0.5
                x = i - sizeInPix/2.0 +0.5
                r = numpy.sqrt(x**2.0+y**2.0)
                f = r/fc
                if f < 1:
                    if r < 0.1:
                        total[i,j] += 1.0/nLambdaSteps
                    else:
                        total[i,j] += (TelOTF(f, M2/M1)*Sinc(numpy.pi*x/sizeInPix)*Sinc(numpy.pi*y/sizeInPix))/nLambdaSteps
                else:
                    total[i,j] += 0.0
    retval = numpy.abs(numpy.fft.fft2(total))
    retval = numpy.fft.fftshift(retval)
    #retval = numpy.real(retval*retval.conj())
    #retval = total
    return retval, total

def Sinc(x):
    if numpy.abs(x) < 1e-4:
        return 1.0
    else:
        return numpy.sin(x)/x

def H1(f, u, v):
    if numpy.abs(1.0-v) < 1e-12:
        e = 1
    else:
        e = -1
    return v**2.0/numpy.pi * numpy.arccos((f/v) * (1+e*(1.0-u**2.0)/(4*f**2.0)))

def H2(f, u):
    a = 2.0*f/(1.0+u)
    b = (1.0-u)/(2.0*f)
    return -1.0*(f/numpy.pi)*(1.0+u)*numpy.sqrt((1-a**2.0)*(1.0-b**2.0))

def G(f, u):
    if f <= (1.0 - u)/2.0:
        return u**2.0
    elif f >= (1.0 + u)/2.0:
        return 0.0
    else:
        return H1(f, u, 1.0) + H1(f, u, u) + H2(f, u)

def TelOTF(f, u):
    retval = (G(f, 1.0) + u**2.0*G(f/u, 1)  - 2*G(f,u))/(1.0-u**2.0)
    return retval


fig = pyplot.figure(0)
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
ax.clear()

image, OTF = generateOTF()

ax.matshow(image)
fig.show()

