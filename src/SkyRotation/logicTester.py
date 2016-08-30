import pyfits
import numpy


derotator = 359.6
azimuth = 139.5

mode = 'STAT'

angleTable = pyfits.getdata('RTC.'+mode+'.ANGLE_TABLE.fits')
HOIM = pyfits.getdata('RTC.'+mode+'.HOIM.fits')

if mode == 'SKY':
    index = numpy.argsort(numpy.abs(angleTable[:,0]-azimuth))[0]
    value = numpy.average(HOIM[:,index])
elif mode == 'STAT':
    maximumAngle = 180.0
    derot = derotator % maximumAngle
    if derot != derotator:
        appliedAngleOffset = maximumAngle
    else:
        appliedAngleOffset = 0.0
    index = numpy.argsort(numpy.abs(angleTable[:,1]-derot))[0]
    minIndex = numpy.min(angleTable[:,1])
    maxIndex = numpy.max(angleTable[:,1])
    if index == maxIndex:
        maxDistance = numpy.abs(derot - angleTable[maxIndex,1])
        minDistance = numpy.abs(derot - maximumAngle - angleTable[minIndex,1])
        if maxDistance > minDistance:
            index = minIndex
            appliedAngleOffset += maximumAngle

    value = numpy.average(HOIM[:,index])

    appliedDerotAngle = value+appliedAngleOffset


print derotator, derot, appliedAngleOffset
print index, value, appliedDerotAngle
