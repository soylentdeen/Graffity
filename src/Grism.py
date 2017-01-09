import scipy
import numpy

class GrismInstrument( object ):
    def __init__(self, Grism = None, f1 = 24.0, f2 = 24.0, CCD = None):
        self.Grism = Grism
        self.f1 = f1
        self.f2 = f2
        self.CCD = CCD

    def setSlitSize(self, imageWidth=2.5, pixelSize=15.0):
        self.deltaX_slit = (imageWidth * self.f1/self.f2)*pixelSize

    def setOpeningAngle(self, Lambda = 2000.0, dLambda = 400.0):
        self.Lambda = Lambda
        self.dLambda = dLambda
        self.Grism.delta = numpy.rad2deg(numpy.arctan(Lambda*self.CCD.nx/self.f2*self.CCD.pitch/1000.0/((self.Grism.n-1.0)*dLambda)))
        self.ResolvingPower = (self.Grism.n - 1.0)*numpy.tan(numpy.deg2rad(self.Grism.delta))*self.f1*1000.0/self.deltaX_slit

    def setSigma(self, order=1):
        self.Grism.sigma = order*self.Lambda/1000.0/((self.Grism.n - 1.0)*numpy.sin(numpy.deg2rad(self.Grism.delta)))

    def optimize(self, npix=2.5, Lambda=2200.0, dLambda=500.0, order=1):
        self.setSlitSize(imageWidth=npix, pixelSize=self.CCD.pitch)
        self.setOpeningAngle(Lambda=Lambda, dLambda=dLambda)
        self.setSigma(order=order)


class Grism( object ):
    def __init__(self, n = 3.41, sigma = 2.5, alpha = 0.0):
        self.n = n
        self.sigma = sigma
        self.alpha = alpha

class CCD( object ):
    def __init__(self, nx = 500, ny = 500, pitch= 18.0):
        self.nx = 2000
        self.ny = ny
        self.pitch = pitch
