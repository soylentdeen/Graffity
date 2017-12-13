import astropy.io.fits as pyfits
from matplotlib import pyplot
import numpy
import glob
import Graffity

datadir = '/home/cdeen/Data/GRAVITY/Raw/reduced/'

fig0 = pyplot.figure(0)
fig0.clear()
ax0 = fig0.add_axes([0.1 ,0.1 , 0.8, 0.8])
fig1 = pyplot.figure(1)
fig1.clear()
ax1 = fig1.add_axes([0.1 ,0.1 , 0.8, 0.8])
fig2 = pyplot.figure(2)
fig2.clear()
ax2 = fig2.add_axes([0.1 ,0.1 , 0.8, 0.8])
fig3 = pyplot.figure(3)
fig3.clear()
ax3 = fig3.add_axes([0.1 ,0.1 , 0.8, 0.8])

#source = "Fringe Tracker"
source = "Science Camera"

size = 20

telescopes = [1, 2, 3, 4]
colors = ['b', 'g', 'r', 'y']

Perfect = Graffity.PSF(sizeInPix=40)
Perfect.generateOTF()

#dark = pyfits.open(datadir+'GRAVI.2017-07-13T07:12:51.032.fits')
#dark_AcqCam = dark["IMAGING_DATA_ACQ"].data[0]
#dark_SC = numpy.median(dark["IMAGING_DATA_SC"].data, axis=0)
#dark_FT = numpy.median(dark["IMAGING_DATA_FT"].data.field("PIX"), axis=0)

AcqCamFlux = {}

FT_table = []
SC_table = []
Acq_FT_table = []
Acq_SC_table = []
for j, f in enumerate(glob.glob(datadir+'*singlescip2vmred.fits')):
    print f
    HDU = pyfits.open(f)
    primary = HDU["PRIMARY"].header
    AcqCam = HDU[17].data[0]
    SC = HDU[10].data.field("TOTALFLUX_FT")
    FT = HDU[10].data.field("TOTALFLUX_SC")


    ax0.clear()
    ax1.clear()
    ax2.clear()
    ax3.clear()
    FT_entry = []
    SC_entry = []
    Acq_FT_entry = []
    Acq_SC_entry = []
    for i in telescopes:
        startX = primary.get("ESO DET1 FRAM%d STRX" %i)
        startY = primary.get("ESO DET1 FRAM%d STRY" %i)
        xcoord_FT = primary.get("ESO ACQ FIBER FT%dX" % i) - startX + (i-1)*250
        ycoord_FT = primary.get("ESO ACQ FIBER FT%dY" % i) - startY
        xcoord_SC = primary.get("ESO ACQ FIBER SC%dX" % i) - startX + (i-1)*250
        ycoord_SC = primary.get("ESO ACQ FIBER SC%dY" % i) - startY
        print ("Telescope %d : Fringe Tracker Flux = %f, Science Camera Flux = %f" % (i, numpy.mean(FT[i-1::4]), numpy.mean(SC[i-1::4])))
        FT_entry.append(numpy.mean(FT[i-1::4]))
        SC_entry.append(numpy.mean(SC[i-1::4]))
        Acq_FT_entry.append(numpy.max(AcqCam[int(ycoord_FT-3):int(ycoord_FT+3), int(xcoord_FT-3):int(xcoord_FT+3)]))
        Acq_SC_entry.append(numpy.max(AcqCam[int(ycoord_SC-3):int(ycoord_SC+3), int(xcoord_SC-3):int(xcoord_SC+3)]))
        #B_entry.append(numpy.sum(AcqCam[int(ycoord_SC-3):int(ycoord_SC+3), int(xcoord_SC-3):int(xcoord_SC+3)]))
        """
        if AcqCam[int(ycoord_FT), int(xcoord_FT)] > 1000.0:
            FT_PostageStamp = AcqCam[int(ycoord_FT-size):int(ycoord_FT+size), int(xcoord_FT-size):int(xcoord_FT+size)]
            FT_Strehl = Perfect.calcStrehl(FT_PostageStamp)
            SC_PostageStamp = AcqCam[int(ycoord_SC-size):int(ycoord_SC+size), int(xcoord_SC-size):int(xcoord_SC+size)]
            SC_Strehl = Perfect.calcStrehl(SC_PostageStamp)
            FTFlux = numpy.sum(AcqCam[int(ycoord_FT-3):int(ycoord_FT+3), int(xcoord_FT-3):int(xcoord_FT+3)])
            SCFlux = numpy.sum(AcqCam[int(ycoord_SC-3):int(ycoord_SC+3), int(xcoord_SC-3):int(xcoord_SC+3)])
            print "Telescope", i, numpy.max(FT_PostageStamp), numpy.max(SC_PostageStamp)
            print ("Telescope %d Fringe Tracking Strehl = %.3f, Science Camera Strehl = %.3f" % (i, FT_Strehl, SC_Strehl))
            ax0.plot(numpy.sum(FT_PostageStamp, axis=0), color = colors[i-1], label="Tel %d SR = %.3f, Flux = %d" % (i, FT_Strehl, FTFlux))
            ax0.set_xlabel("AcqCam FT Object" % FT_Strehl)
            ax1.plot(numpy.sum(SC_PostageStamp, axis=0), color = colors[i-1], label="Tel %d SR = %.3f, Flux = %d" % (i, SC_Strehl, SCFlux))
            ax1.set_xlabel("AcqCam SC Object" % SC_Strehl)
        #"""
    if not (j in [4, 7, 8]):
        SC_table.append(numpy.array(SC_entry))
        FT_table.append(numpy.array(FT_entry))
        Acq_FT_table.append(numpy.array(Acq_FT_entry))
        Acq_SC_table.append(numpy.array(Acq_SC_entry))

SC_table = numpy.array(SC_table)
Acq_SC_table = numpy.array(Acq_SC_table)
FT_table = numpy.array(FT_table)
Acq_FT_table = numpy.array(Acq_FT_table)

x_FT = numpy.linalg.pinv(Acq_FT_table)
x_SC = numpy.linalg.pinv(Acq_SC_table)

ans_FT = x_FT.dot(FT_table)
ans_SC = x_SC.dot(SC_table)

if source == 'Fringe Tracker':
    ax0.matshow(Acq_FT_table)
    cov = ax1.matshow(ans_FT)
    fig1.colorbar(cov)
    ax2.matshow(FT_table)
    ax3.matshow(Acq_FT_table.dot(ans_FT))
elif source == 'Science Camera':
    ax0.matshow(Acq_SC_table)
    cov = ax1.matshow(ans_SC)
    fig1.colorbar(cov)
    ax2.matshow(SC_table)
    ax3.matshow(Acq_SC_table.dot(ans_SC))

ax0.set_title("AcqCamera %s - Max Pixel" % source)
ax0.set_xlabel("TEL channel")
ax0.set_ylabel("Trial")
ax1.set_title("%s Covariance Matrix" % source)
ax2.set_title("%s Pipeline Flux" % source)
ax2.set_xlabel("TOTALFLUX channel")
ax2.set_ylabel("Trial")
ax3.set_title("%s Fit" % source)
ax3.set_xlabel("TOTALFLUX channel")
ax3.set_ylabel("Trial")
fig0.show()
fig1.show()
fig2.show()
fig3.show()
print "Science Camera: ", ans_SC.diagonal()
print "Fringe Tracker: ", ans_FT.diagonal()

fig0.savefig("AcqCam_%s_MaxPix.png" % source.replace(' ', '_'))
fig1.savefig("%s_CovMat.png" % source.replace(' ', '_'))
fig2.savefig("%s_Pipeline_Flux.png" % source.replace(' ', '_'))
fig3.savefig("%s_fit.png" % source.replace(' ', '_'))



"""
    #ax2.plot(numpy.sum(FT, axis=0))
    #ax2.set_xlabel("Fringe Tracker")
    #ax3.plot(numpy.sum(SC, axis=1))
    #ax3.set_xlabel("Science Object")
    #print asddf
    ax1.legend()
    ax0.legend()
    fig0.show()
    fig1.show()
    #fig2.show()
    #fig3.show()
    decision = raw_input("Save? ")
    if decision== 'Y':
        baseName = raw_input("Enter base filename: ")
        fig0.savefig(baseName+"_FT_AcqCam.png")
        fig1.savefig(baseName+"_SC_AcqCam.png")
        fig2.savefig(baseName+"_FT_Trace.png")
        fig3.savefig(baseName+"_SC_Trace.png")
    #"""
