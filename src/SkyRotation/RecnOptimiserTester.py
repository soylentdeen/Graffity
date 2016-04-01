import pyfits
import numpy

STAT_HOIM = []
STAT_ANGSENS = []
STAT_XSENS = []
STAT_YSENS = []
STAT_ANGLE_TABLE = []
for derotAngle in numpy.arange(180.0):
    IM = numpy.ones(8160)*(-1.0*derotAngle)
    STAT_HOIM.append(IM)
    ANGSENS = (numpy.ones(8160)*derotAngle+0.1)*(-1.0)
    STAT_ANGSENS.append(ANGSENS)
    XSENS = (numpy.ones(8160)*derotAngle+0.2)*(-1.0)
    STAT_XSENS.append(XSENS)
    YSENS = (numpy.ones(8160)*derotAngle+0.3)*(-1.0)
    STAT_YSENS.append(YSENS)
    ANGLE_TABLE = numpy.array([139.5, derotAngle])
    STAT_ANGLE_TABLE.append(ANGLE_TABLE)

STAT_HOIM = numpy.array(STAT_HOIM, dtype=numpy.float32)
STAT_ANGSENS = numpy.array(STAT_ANGSENS, dtype=numpy.float32)
STAT_YSENS = numpy.array(STAT_YSENS, dtype=numpy.float32)
STAT_XSENS = numpy.array(STAT_XSENS, dtype=numpy.float32)
STAT_ANGLE_TABLE = numpy.array(STAT_ANGLE_TABLE, dtype=numpy.float32)

hdu = pyfits.PrimaryHDU(data=STAT_HOIM.T)
hdu.writeto('RTC.STAT.HOIM.fits', clobber=True)
hdu = pyfits.PrimaryHDU(data=STAT_ANGSENS.T)
hdu.writeto('RTC.STAT.ANGSENS.fits', clobber=True)
hdu = pyfits.PrimaryHDU(data=STAT_XSENS.T)
hdu.writeto('RTC.STAT.XSENS.fits', clobber=True)
hdu = pyfits.PrimaryHDU(data=STAT_YSENS.T)
hdu.writeto('RTC.STAT.YSENS.fits', clobber=True)
hdu = pyfits.PrimaryHDU(data=STAT_ANGLE_TABLE)
hdu.writeto('RTC.STAT.ANGLE_TABLE.fits', clobber=True)


SKY_HOIM = []
SKY_ANGSENS = []
SKY_XSENS = []
SKY_YSENS = []
SKY_ANGLE_TABLE = []
for skyAngle in numpy.arange(0.0, 360.0, 2.0):
    IM = numpy.ones(8160)*(skyAngle)
    SKY_HOIM.append(IM)
    ANGSENS = numpy.ones(8160)*skyAngle+0.1
    SKY_ANGSENS.append(ANGSENS)
    XSENS = numpy.ones(8160)*skyAngle +0.2
    SKY_XSENS.append(XSENS)
    YSENS = numpy.ones(8160)*skyAngle +0.3
    SKY_YSENS.append(YSENS)
    ANGLE_TABLE = numpy.array([skyAngle, 0.0])
    SKY_ANGLE_TABLE.append(ANGLE_TABLE)

SKY_HOIM = numpy.array(SKY_HOIM, dtype=numpy.float32)
SKY_ANGSENS = numpy.array(SKY_ANGSENS, dtype=numpy.float32)
SKY_YSENS = numpy.array(SKY_YSENS, dtype=numpy.float32)
SKY_XSENS = numpy.array(SKY_XSENS, dtype=numpy.float32)
SKY_ANGLE_TABLE = numpy.array(SKY_ANGLE_TABLE, dtype=numpy.float32)

hdu = pyfits.PrimaryHDU(data=SKY_HOIM.T)
hdu.writeto('RTC.SKY.HOIM.fits', clobber=True)
hdu = pyfits.PrimaryHDU(data=SKY_ANGSENS.T)
hdu.writeto('RTC.SKY.ANGSENS.fits', clobber=True)
hdu = pyfits.PrimaryHDU(data=SKY_XSENS.T)
hdu.writeto('RTC.SKY.XSENS.fits', clobber=True)
hdu = pyfits.PrimaryHDU(data=SKY_YSENS.T)
hdu.writeto('RTC.SKY.YSENS.fits', clobber=True)
hdu = pyfits.PrimaryHDU(data=SKY_ANGLE_TABLE)
hdu.writeto('RTC.SKY.ANGLE_TABLE.fits', clobber=True)

