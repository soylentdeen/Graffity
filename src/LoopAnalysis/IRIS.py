import scipy
import numpy
from matplotlib import pyplot
import Graffity
import glob
import time
import astropy.io.fits as pyfits
from scipy import ndimage

#IRIS_Data = pyfits.getdata('/home/cdeen/Data/CIAO/JanComm/IRIS20170112/hd_22468_low_DIT.fits') - \
#                  numpy.mean(Background, axis=0)
#IRIS_Data = pyfits.getdata('/home/cdeen/Data/CIAO/JanComm/IRIS20170112/v368_pup_low_DIT.fits') - \
#                  numpy.mean(Background, axis=0)
#IRIS_Data = pyfits.getdata('/home/cdeen/Data/CIAO/JanComm/IRIS20170112/wds_j03575_DIT.fits') - \
#                  numpy.mean(Background, axis=0)

def calcStrehl(base='', sens='', airy='', blah=False):
    if blah:
        base_dir = '/home/cdeen/Data/CIAO/JanComm/IRIS20170108/'
        Background = pyfits.getdata(base_dir+'HD41_BACKGROUND_1_DIT.fits')
        if sens == 'high':
            IRIS_Data = pyfits.getdata(base_dir+'HD41_AVC_ON_1_DIT.fits')-numpy.mean(Background, axis=0)
            IRIS_Header = pyfits.getheader(base_dir+'HD41_AVC_ON_1_DIT.fits')
        else:
            IRIS_Data = pyfits.getdata(base_dir+'HD41_faint_AVC_ON_1_DIT.fits')-numpy.mean(Background, axis=0)
            IRIS_Header = pyfits.getheader(base_dir+'HD41_faint_AVC_ON_1_DIT.fits')

        print IRIS_Header.get("DATE-OBS")
    else:
        base_dir = '/home/cdeen/Data/CIAO/JanComm/IRIS20170112/'
        Background=pyfits.getdata(base_dir+'iris_background_5ms_DIT.fits')
        IRIS_Data = pyfits.getdata(base_dir + base + '_' + sens + '_DIT.fits') - numpy.mean(Background, axis=0)

    UT1 = IRIS_Data[:, :32, :]
    UT2 = IRIS_Data[:, 32:64, :]
    UT3 = IRIS_Data[:, 64:96, :]
    UT4 = IRIS_Data[:, 96:, :]

    #UT1 = numpy.mean(UT1, axis=0)
    #UT2 = numpy.mean(UT2, axis=0)
    #UT3 = numpy.mean(UT3, axis=0)
    #UT4 = numpy.mean(UT4, axis=0)


    SR = {}
    centroid = {}
    for UT, stack in zip([1, 2, 3, 4], [UT1, UT2, UT3, UT4]):
        SR[UT] = []
        centroid[UT] = []
        for frame in stack:
            frame = frame / numpy.max(frame)
            frame[frame < 0] = 0.0
            SR[UT].append(numpy.max(frame/numpy.sum(frame)*numpy.sum(airy)))
            centroid[UT].append(ndimage.measurements.center_of_mass(frame))

        SR[UT] = numpy.array(SR[UT])
        centroid[UT] = numpy.std(numpy.array(centroid[UT]), axis=0)*31.5
    
    return SR, centroid

fig0 = pyplot.figure(0)
fig0.clear()
ax0 = fig0.add_axes([0.1, 0.1, 0.8, 0.8])
fig1 = pyplot.figure(1)
fig1.clear()
ax1 = fig1.add_axes([0.1, 0.1, 0.8, 0.8])

PSF = Graffity.PSF(sizeInPix=32, lam=2.2, pscale = 31.5)
PSF.generateOTF()
airy = PSF.getPSF()

#stars = ['v368_pup', 'hd_22468', 'wds_j03575']
stars = ['v368_pup', 'wds_j03575']
sensitivities = ['low', 'medium', 'high']
#separations = [11.8, 6.7, 11.1]
separations = [11.8, 11.1]

SR = {}
TTJit = {}
for star, sep in zip(stars, separations):
    SR[star] = {}
    SR[star]["separation"] = sep
    TTJit[star] = {}
    TTJit[star] = {}
    for sens in sensitivities:
        blah = calcStrehl(base=star, sens=sens, airy=airy)
        SR[star][sens] = blah[0]
        TTJit[star][sens] = blah[1]
        print("%s  %s" % (star, sens))
        print("CIAO 1: %.3f" % (numpy.mean(SR[star][sens][1])))
        print("CIAO 2: %.3f" % (numpy.mean(SR[star][sens][2])))
        print("CIAO 3: %.3f" % (numpy.mean(SR[star][sens][3])))
        print("CIAO 4: %.3f" % (numpy.mean(SR[star][sens][4])))

low, ttj_low = calcStrehl(sens='low', airy=airy, blah = True)
high, ttj_high = calcStrehl(sens='high', airy=airy, blah = True)

SR["HD41"] = {}
SR["HD41"]["separation"] = 5.58
SR["HD41"]["medium"] = low
SR["HD41"]["high"] = high
SR["HD41"]["low"] = low

TTJit["HD41"] = {}
TTJit["HD41"]["separation"] = 5.58
TTJit["HD41"]["medium"] = ttj_low
TTJit["HD41"]["high"] = ttj_high
TTJit["HD41"]["low"] = ttj_low

separations = []
strehl_1 = []
strehl_2 = []
strehl_3 = []
strehl_4 = []
ttx_1 = []
tty_1 = []
ttx_2 = []
tty_2 = []
ttx_3 = []
tty_3 = []
ttx_4 = []
tty_4 = []
sens = []
for star in SR.keys():
    for s in sensitivities:
        separations.append(SR[star]["separation"])
        strehl_1.append(numpy.mean(SR[star][s][1]))
        strehl_2.append(numpy.mean(SR[star][s][2]))
        strehl_3.append(numpy.mean(SR[star][s][3]))
        strehl_4.append(numpy.mean(SR[star][s][4]))
        ttx_1.append(TTJit[star][s][1][0])
        tty_1.append(TTJit[star][s][1][1])
        ttx_2.append(TTJit[star][s][2][0])
        tty_2.append(TTJit[star][s][2][1])
        ttx_3.append(TTJit[star][s][3][0])
        tty_3.append(TTJit[star][s][3][1])
        ttx_4.append(TTJit[star][s][4][0])
        tty_4.append(TTJit[star][s][4][1])

ax0.scatter(separations, ttx_1, color = 'b', marker='x', s=40.0)
ax0.scatter(numpy.array(separations)+0.1, tty_1, color = 'b', marker='+', s=40.0)
ax0.scatter(separations, ttx_2, color = 'g', marker='x', s=40.0)
ax0.scatter(numpy.array(separations)+0.1, tty_2, color = 'g', marker='+', s=40.0)
ax0.scatter(separations, ttx_3, color = 'r', marker='x', s=40.0)
ax0.scatter(numpy.array(separations)+0.1, tty_3, color = 'r', marker='+', s=40.0)
ax0.scatter(separations, ttx_4, color = 'm', marker='x', s=40.0)
ax0.scatter(numpy.array(separations)+0.1, tty_4, color = 'm', marker='+', s=40.0)

ax0.set_title("UT1 - Tip/Tilt Jitter vs Separation")
ax0.set_ylabel("Tip/Tilt Jitter (mas)")

#ax0.scatter(separations, strehl_1, color='b')
#ax0.scatter(separations, strehl_2, color='g')
#ax0.scatter(separations, strehl_3, color='r')
#ax0.scatter(separations, strehl_4, color='m')
#ax0.set_title("UT1 - Strehl vs Separation")
ax0.set_xlabel("Separation - arcsec")
#ax0.set_ylabel("Strehl Ratio (IRIS)")

#ax0.matshow(airy)
#ax1.matshow(numpy.mean(UT4, axis=0))

fig0.show()
fig0.savefig("TTJ_v_Separation.png")
#fig0.savefig("SR_v_Separation.png")
#fig1.show()
