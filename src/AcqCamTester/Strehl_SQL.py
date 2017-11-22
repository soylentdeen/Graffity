import astropy.io.fits as pyfits
from matplotlib import pyplot
import numpy
import scipy
from scipy.stats import kde
import glob
import sys

sys.path.append('../')
import CIAO_DatabaseTools
import Graffity
import astropy.time as aptime

def findCorrelations(GravObs=[], ax=None):
    colors = ['b', 'g', 'r', 'y']
    colormaps = [pyplot.cm.Blues_r, pyplot.cm.Greens_r, pyplot.cm.Reds_r,
            pyplot.cm.Oranges_r]
    #if ax != None:
    #    ax[0].clear()
    #    ax[1].clear()
    #    ax[2].clear()
    linex = numpy.array([numpy.ones(50), numpy.linspace(0, 0.1)])
    difference = {}
    strehl = {}
    cov = {}
    CIAO_SR = {}
    ellipsicity = {}
    baselines = {0:[4,3], 1:[4, 2], 2:[4, 1], 3:[3, 2], 4:[3, 1], 5:[2, 1]}

    for i,j in zip(range(4), [3, 2, 1,0]):
        A = []
        A_FT = []   #Fringe Tracker
        A_SC = [] #Science Object
        x = []
        x_FT = []   #Fringe Tracker
        x_SC = []   #Science Object
        strehl[i] = []
        on_error = []
        off_error = []
        Pupil = []
        CIAO_SR[i] = []
        ellipsicity[i] = []
        for Grav in GravObs:
            DimTripwire = True
            BrightTripwire = True
            #if hasattr(Grav.DualSciP2VM, "AcqCamDat"):
            derot = Grav.getDerotatorPositions()
            if Grav.AcqCamData != None:
                Acq = Grav.AcqCamData
                strehlAccumulator = []
                if Acq.CIAO_Data != None:
                    CIAO_SR[i].append(Acq.CIAO_Data[j]["Strehl"])
                fluxAccumulator = []
                Grav.computeOPDPeriodograms()
                peaks = Grav.findVibrationPeaks()
                for SR, FT_flux, SC_flux, dx, dy, sc, pupil in zip(Acq.newStrehl.data[i],
                        Acq.TOTALFLUX_FT.data[i], Acq.TOTALFLUX_SC.data[i], Acq.newSC_FIBER_DX.data[i],
                        Acq.newSC_FIBER_DY.data[i], Acq.newSCALE.data[i],
                        Acq.SC_Pos.data[i]):
                    if (not(numpy.isnan(SR)) and (SR > 0) and (SR < 0.5)):
                        if (FT_flux > 150000):
                            strehl[i].append(SR)
                            #ellipsicity[i].append(ellipse)

                            strehlAccumulator.append(SR)
                            fluxAccumulator.append(SC_flux)
                            #A.append([1.0, SR/numpy.max(Acq.newStrehl.data[i])])
                            #x.append(flux / Acq.TOTALFLUX_FT.median[i])
                            #A.append([1.0, SR])
                            A.append([1.0, SR, 10.0**(-2.5*Acq.FTMag)])
                            A.append([1.0, SR, 10.0**(-2.5*Acq.SCMag)])
                            A_FT.append([1.0, SR])
                            A_SC.append([1.0, SR])
                            #A.append([1.0, SR, 10.0**(-2.5*Acq.FTMag)*Acq.AcqDit,
                            #    Acq.CIAO_Data[j]["Strehl"],
                            #    Acq.CIAO_Data[j]["Seeing"],
                            #    Acq.CIAO_Data[j]["Tau0"]])
                            #    Acq.CIAO_Data[j]["WindSp"]])
                            #A.append([1.0, SR, (dx**2.0+dy**2.0)**0.5])
                            #A.append([1.0, SR, dx, dy])
                            #A.append([1.0, SR, dx, dy, sc])
                            x_FT.append(FT_flux)
                            x_SC.append(SC_flux)
                            x.append(FT_flux)
                            x.append(SC_flux)
                            Pupil.append(pupil)
                            on_error.append((dx**2.0+dy**2.0)**0.5)
                            #if BrightTripwire:
                            #    print("Good Filename = %s" % Grav.DualSciP2VM.filename)
                            #    BrightTripwire = False
                            """
                            if ax != None:
                                ax[2][i].scatter([SC_flux], [dx**2.0 + dy**2.0],
                                        c=colors[i], marker='.')
                            #"""
                        else:
                            if not(numpy.isnan(dx) or numpy.isnan(dy)):
                                off_error.append((dx**2.0+dy**2.0)**0.5)
                            #if DimTripwire:
                            #    print("Bad Filename = %s" % Grav.DualSciP2VM.filename)
                            #    DimTripwire = False
                            #print("dx = %.3f, dy = %.3f" %
                            #        (dx, dy))
                if ax!= None:
                    #ax[2][i].scatter(numpy.median(fluxAccumulator),
                    #        [Grav.startTime.mjd - 57970.0], c=colors[i])
                    #ax[2][i].scatter([Grav.startTime.mjd - 57970.0],
                    #        [numpy.median(fluxAccumulator)], c=colors[i])
                    ratio = numpy.median(numpy.array(fluxAccumulator)/numpy.array(strehlAccumulator))
                    #ax[3].scatter([derot[i]-derot[i+4]], [ratio],
                    ax[3].scatter([derot[i]],
                            [numpy.median(numpy.array(fluxAccumulator))],
                            c=colors[i])
                affected = []
                for bl in baselines.keys():
                    if j+1 in baselines[bl]:
                        affected.append(bl)

                freqs = []
                for bl in affected:
                    freqs.append(peaks[bl]['freqs'])

                uniqueFreqs = numpy.unique(numpy.append(numpy.append(freqs[0],
                    freqs[1]), freqs[2]))
                CommonFreqs = []
                for uniq in uniqueFreqs:
                    if ((uniq in freqs[0]) and (uniq in freqs[1]) and (uniq in
                        freqs[2])):
                        CommonFreqs.append(uniq)
                
                for comm in CommonFreqs:
                    for bl in affected:
                        ax[4].scatter([peaks[bl]['power'][peaks[bl]['freqs']==comm]],[numpy.mean(strehlAccumulator)],color=colors[i])



                                



        A = numpy.array(A, dtype=float)
        x = numpy.array(x)
        A_FT = numpy.array(A_FT, dtype=float)
        x_FT = numpy.array(x_FT)
        A_SC = numpy.array(A_SC, dtype=float)
        x_SC = numpy.array(x_SC)
        B = numpy.linalg.pinv(A)
        fit = B.dot(x)
        B_FT = numpy.linalg.pinv(A_FT)
        fit_FT = B_FT.dot(x_FT)
        B_SC = numpy.linalg.pinv(A_SC)
        fit_SC = B_SC.dot(x_SC)
        print fit
        print fit_FT, A_FT.shape
        print fit_SC, A_SC.shape
        print ("On error = %.3f" % numpy.mean(numpy.array(on_error)))
        print ("Off error = %.3f" % numpy.mean(numpy.array(off_error)))
        #ax[0].scatter(A_FT[:,1], x_FT, c=colors[i])
        counts, yedges, xedges, image = ax[0][i].hist2d(A_FT[:,1], x_FT, bins=50,
                cmap=colormaps[i])
        #nbins = 20
        #xi, yi = numpy.mgrid[xedges.min():xedges.max():nbins*1j,
        #    yedges.min():yedges.max():nbins*1j]
        #k = kde.gaussian_kde(numpy.array([A_FT[:,1], x_FT]))
        #zi = k(numpy.vstack([xi.flatten(), yi.flatten()]))
        #ax[0].pcolormesh(yi, xi, zi.reshape(xi.shape), shading='gouraud',
        #        cmap=colormaps[i], alpha=0.3)
        #ax[0].contour(counts, extent=[xedges.min(), xedges.max(), yedges.min(),
        #    yedges.max()],cmap=colormaps[i], alpha=0.5,
        #    levels=[counts.max()/5.0])
        #ax[1].scatter(A_SC[:,1], x_SC, c=colors[i])
        counts, yedges, xedges, image = ax[1][i].hist2d(A_SC[:,1], x_SC, bins=50,
                cmap=colormaps[i])
        #junk = numpy.hstack((linex.T, numpy.ones([linex.shape[1],
        #    1])*10.0**(-2.5*Acq.FTMag)))
        CIAO_SR[i] = numpy.array(CIAO_SR[i], dtype=float)
        strehl[i] = numpy.array(strehl[i])
        ellipsicity[i] = numpy.array(ellipsicity[i])
        for n in range(4):
            ax[0][n].plot(linex[1,:], fit_FT.dot(linex), color = colors[i])
            ax[1][n].plot(linex[1,:], fit_SC.dot(linex), color = colors[i])
        #ax[2][i].hist(CIAO_SR[i], bins=numpy.linspace(0,1.0, num=10),
        #            color=colors[i])
        #ax[4].scatter(ellipsicity[i], x_FT, color=colors[i])
        ax[2][i].hist(numpy.array(strehl[i], dtype=float),
                bins=numpy.linspace(0,0.35),
                    color=colors[i])
        #print("CIAO Strehl")
        #print("Mean = %.3f Stdev = %.3f Skew = %.3f" % (numpy.mean(CIAO_SR[i]),
        #    numpy.std(CIAO_SR[i]), scipy.stats.skew(CIAO_SR[i])))
        print("AcqCam Strehl")
        print("Mean = %.3f Stdev = %.3f Skew = %.3f" % (numpy.mean(strehl[i]),
            numpy.std(strehl[i]), scipy.stats.skew(strehl[i])))

        #n = A.shape[0]
        #avg = numpy.mean(A, axis=0)
        #cov[i] = 1.0/n * (A - avg).T.dot(A-avg)

    for i in range(4):
        #print("KS Test for CIAO and AcqCam SR for Tel%d: %.3f" % (i+1,
        #    scipy.stats.ks_2samp(CIAO_SR[i], strehl[i]).pvalue))
        for j in range(i):
            if i != j:
                #print("KS Test for CIAO SR for Tel%d and Tel%d: %.3f" % (i+1, j+1,
                #    scipy.stats.ks_2samp(CIAO_SR[i], CIAO_SR[j]).pvalue))
                print("KS Test for AcqCam SR for Tel%d and Tel%d: %.3f" % (i+1, j+1,
                    scipy.stats.ks_2samp(strehl[i], strehl[j]).pvalue))

    #ax[0].set_xlabel("H-band AcqCam Strehl")
    #ax[0].set_ylabel("FT Flux")
    #ax[0].set_title("IRS16C, August")
    #ax[0].set_title("IRS16C, July")
    #ax[1].set_xlabel("H-band AcqCam Strehl")
    #ax[1].set_ylabel("SC Flux")
    #ax[1].set_title("S2, August")
    #ax[1].set_title("S2, July")
    #if ax != None:
    #    ax[1].clear()
    #    ax[1].plot(strehl[0], c=colors[0])
    #    ax[1].plot(strehl[1], c=colors[1])
    #    ax[1].plot(strehl[2], c=colors[2])
    #    ax[1].plot(strehl[3], c=colors[3])

    return "cov"

def getCIAO_DataLogger(Grav, DataLoggers):
    retval = {}
    for i in [1, 2, 3, 4]:
        retval[i-1] = {}
        times = numpy.array(DataLoggers[i][:,-4], dtype=float)
        closest = numpy.argsort(numpy.abs(times - float(Grav[-2])))[0]
        if (numpy.abs(times[closest] - float(Grav[-2])) >
                (1.0/(24.0*3600/30.0))):
            good = False
            return None
        else:
            DL = Graffity.DataLogger(CIAO_DataLoggers[i][closest,-3])
            DL.loadData()
            DL.getRefSlopeZernikes()
            retval[i-1]["Tau0"] = CIAO_DataLoggers[i][closest,0]
            retval[i-1]["Strehl"] = CIAO_DataLoggers[i][closest,1]
            retval[i-1]["Seeing"] = CIAO_DataLoggers[i][closest,2]
            retval[i-1]["WindSp"] = CIAO_DataLoggers[i][closest,3]
            retval[i-1]["RefSlopes"] = DL.rotatedZernikes
            del(DL)
            
    return retval


fig1 = pyplot.figure(1)
fig1.clear()
ax11 = fig1.add_axes([0.1, 0.1, 0.4, 0.4])
ax12 = fig1.add_axes([0.1, 0.5, 0.4, 0.4], sharex=ax11, sharey=ax11)
ax13 = fig1.add_axes([0.5, 0.1, 0.4, 0.4], sharey=ax11, sharex=ax11)
ax14 = fig1.add_axes([0.5, 0.5, 0.4, 0.4], sharex=ax13, sharey=ax12)
ax11.clear()
ax12.clear()
ax13.clear()
ax14.clear()
fig2 = pyplot.figure(2)
fig2.clear()
ax21 = fig2.add_axes([0.1, 0.1, 0.4, 0.4])
ax22 = fig2.add_axes([0.1, 0.5, 0.4, 0.4], sharex=ax21, sharey=ax21)
ax23 = fig2.add_axes([0.5, 0.1, 0.4, 0.4], sharey=ax21, sharex=ax21)
ax24 = fig2.add_axes([0.5, 0.5, 0.4, 0.4], sharex=ax23, sharey=ax22)
ax21.clear()
ax22.clear()
ax23.clear()
ax24.clear()
fig3 = pyplot.figure(3)
fig3.clear()
ax31 = fig3.add_axes([0.1, 0.1, 0.4, 0.4])
ax32 = fig3.add_axes([0.1, 0.5, 0.4, 0.4])
ax33 = fig3.add_axes([0.5, 0.1, 0.4, 0.4])
ax34 = fig3.add_axes([0.5, 0.5, 0.4, 0.4])
ax31.clear()
ax32.clear()
ax33.clear()
ax34.clear()
fig4 = pyplot.figure(4)
fig4.clear()
ax4 = fig4.add_axes([0.1, 0.1, 0.8, 0.8])
fig5 = pyplot.figure(5)
fig5.clear()
ax5 = fig5.add_axes([0.1, 0.1, 0.8, 0.8])

AcqFiles = []

GDB = CIAO_DatabaseTools.GRAVITY_Database()
CDB = CIAO_DatabaseTools.CIAO_Database()
CIAO_DataLoggers = CDB.query(keywords = ["TAU0", "STREHL", "SEEING", "WINDSP"])

startTime = '2017-01-01 00:00:00'
stopTime = '2018-09-01 00:00:00'

GravityVals = GDB.query(keywords = ['FTOBJ_NAME', 'SOBJ_NAME', 'FTMAG',
        'SOBJMAG', 'AO_SYSTEM', 'DEROT1', 'DEROT2', 'DEROT3', 'DEROT4',
        'DEROT5', 'DEROT6', 'DEROT7', 'DEROT8'], timeOfDay='NIGHT', startTime=startTime,
        endTime=stopTime)

GravData = []
for Grav in GravityVals:
    CIAO_Data = getCIAO_DataLogger(Grav, CIAO_DataLoggers)
    #if ((CIAO_Data != None)):
    #if ((CIAO_Data != None) and (Grav[1] == 'S2') and (Grav[0] == 'IRS16C')):
    #if ((Grav[4] == 'MACAO') and (CIAO_Data==None) and (Grav[0] == 'SS433')):
    #if ((Grav[0] == '3C_273')):
    if ((Grav[4] == 'CIAO') and (Grav[1]=='S2') and (Grav[0] == 'IRS16C')):
        GravData.append(Graffity.GRAVITY_Data(fileBase=Grav[-1],
            CIAO_Data=CIAO_Data, processAcqCamData=False))
        print Grav[-1]

        
covar = findCorrelations(GravData, ax=[[ax11, ax12, ax13, ax14], [ax21, ax22,
    ax23, ax24], [ax31, ax32, ax33, ax34], ax4, ax5])

#Acq.binData()
#Acq.plot(axes=axes)
#Acq.findCorrelations(ax1)

ax11.set_xbound(0.0, 0.15)
ax11.set_ybound(ax11.dataLim.get_points()[0][1], ax11.dataLim.get_points()[1][1])
ax11.set_xlabel("Strehl Ratio")
ax11.set_ylabel("Fringe Tracking Flux")
pyplot.setp(ax12.get_xticklabels(), visible=False)
pyplot.setp(ax14.get_xticklabels(), visible=False)
pyplot.setp(ax13.get_yticklabels(), visible=False)
pyplot.setp(ax14.get_yticklabels(), visible=False)
ax21.set_xbound(0.0, 0.15)
ax21.set_ybound(ax21.dataLim.get_points()[0][1], ax21.dataLim.get_points()[1][1])
ax21.set_xlabel("Strehl Ratio")
ax21.set_ylabel("Science Object Flux")
pyplot.setp(ax22.get_xticklabels(), visible=False)
pyplot.setp(ax24.get_xticklabels(), visible=False)
pyplot.setp(ax23.get_yticklabels(), visible=False)
pyplot.setp(ax24.get_yticklabels(), visible=False)

ax31.set_xlabel("Strehl Ratio")
ax31.set_ylabel("N")

ax4.set_xlabel("CIAO Strehl Ratio")
ax4.set_ylabel("AcqCam Strehl Ratio")
fig4.suptitle("AcqCam vs CIAO Strehl Ratio")
#ax31.set_xbound(0, 2.0)
#ax32.set_xbound(0, 2.0)
#ax33.set_xbound(0, 2.0)
#ax34.set_xbound(0, 2.0)

fig1.suptitle("Fringe Tracking Object = IRS16C")
fig2.suptitle("Science Object = S2")
fig1.show()
fig2.show()
fig3.show()
fig4.show()
fig5.show()
fig1.savefig("FT_Flux.png")
fig2.savefig("SC_Flux.png")
fig3.suptitle("CIAO Strehl Ratios")
#fig3.savefig("Time.png")
#fig3.suptitle("Fiber Positioning Error")
#fig3.savefig("CIAOStrehl.png")
fig3.savefig("CIAOStrehl.png")
fig4.savefig("CIAO_vs_AcqCam_Strehl.png")
#fig1.savefig("FT_August.png")
#fig2.savefig("SC_August.png")
#fig3.show()
#fig4.show()
#fig5.show()
#fig6.show()
