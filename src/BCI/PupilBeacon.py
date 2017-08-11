import Graffity
import scipy
from matplotlib import pyplot
from matplotlib.backends.backend_pdf import PdfPages

fig1 = pyplot.figure(1)
ax11 = fig1.add_axes([0.1, 0.1, 0.4, 0.4])
ax12 = fig1.add_axes([0.1, 0.5, 0.4, 0.4])
ax13 = fig1.add_axes([0.5, 0.1, 0.4, 0.4])
ax14 = fig1.add_axes([0.5, 0.5, 0.4, 0.4])
fig2 = pyplot.figure(2)
ax21 = fig2.add_axes([0.1, 0.1, 0.4, 0.4])
ax22 = fig2.add_axes([0.1, 0.5, 0.4, 0.4])
ax23 = fig2.add_axes([0.5, 0.1, 0.4, 0.4])
ax24 = fig2.add_axes([0.5, 0.5, 0.4, 0.4])
ax21.clear()
ax22.clear()
ax23.clear()
ax24.clear()

datadir= '/home/cdeen/Data/GRAVITY/2017-08/'

for prefix in ['NoAcq_NoBTK', 'Acq_NoBTK', 'Acq_BTK']:
    pp = PdfPages(prefix+'_PSDs.pdf')
    for TT in [50, 100, 200, 400, 800, 1600]:
        df = datadir+'TT'+str(TT)+'_'+prefix+'_20170803.fits'

        title = 'TT'+str(TT)+'_'+prefix
        print title

        LBD = Graffity.LaserBeaconData(df)
        LBD.computePSDs()
        ax11.clear()
        ax12.clear()
        ax13.clear()
        ax14.clear()
        ax21.clear()
        ax22.clear()
        ax23.clear()
        ax24.clear()

        ax11.set_yscale('log')
        ax11.set_xscale('linear')
        ax12.set_yscale('log')
        ax12.set_xscale('linear')
        ax13.set_yscale('log')
        ax13.set_xscale('linear')
        ax14.set_yscale('log')
        ax14.set_xscale('linear')
        ax11.text(400., 1e-5, 'UT1')
        ax12.text(400., 1e-5, 'UT2')
        ax13.text(400., 1e-5, 'UT3')
        ax14.text(400., 1e-5, 'UT4')
        ax11.plot(LBD.Freq, LBD.PSDPower[1])
        ax12.plot(LBD.Freq, LBD.PSDPower[2])
        ax13.plot(LBD.Freq, LBD.PSDPower[3])
        ax14.plot(LBD.Freq, LBD.PSDPower[4])
        ax11.set_ybound(1e-3, 1e-11)
        ax11.set_xbound(0.0, 500)
        ax12.set_ybound(1e-3, 1e-11)
        ax12.set_xbound(0.0, 500)
        ax13.set_ybound(1e-3, 1e-11)
        ax13.set_xbound(0.0, 500)
        ax14.set_ybound(1e-3, 1e-11)
        ax14.set_xbound(0.0, 500)
        ax12.set_xticklabels([])
        ax13.yaxis.tick_right()
        ax14.set_xticklabels([])
        ax14.yaxis.tick_right()

        ax21.set_xscale('linear')
        ax21.set_yscale('log')
        ax22.set_xscale('linear')
        ax22.set_yscale('log')
        ax23.set_xscale('linear')
        ax23.set_yscale('log')
        ax24.set_xscale('linear')
        ax24.set_yscale('log')
        ax21.text(400., 1e-5, 'UT1')
        ax22.text(400., 1e-5, 'UT2')
        ax23.text(400., 1e-5, 'UT3')
        ax24.text(400., 1e-5, 'UT4')
        ax21.plot(LBD.Freq, LBD.PiezoPower[1][:,0:2])
        ax22.plot(LBD.Freq, LBD.PiezoPower[2][:,0:2])
        ax23.plot(LBD.Freq, LBD.PiezoPower[3][:,0:2])
        ax24.plot(LBD.Freq, LBD.PiezoPower[4][:,0:2])
        ax21.set_ybound(1e-3, 1e-11)
        ax21.set_xbound(0.0, 500)
        ax22.set_ybound(1e-3, 1e-11)
        ax22.set_xbound(0.0, 500)
        ax23.set_ybound(1e-3, 1e-11)
        ax23.set_xbound(0.0, 500)
        ax24.set_ybound(1e-3, 1e-11)
        ax24.set_xbound(0.0, 500)

        ax22.set_xticklabels([])
        ax23.yaxis.tick_right()
        ax24.set_xticklabels([])
        ax24.yaxis.tick_right()

        fig1.suptitle(title+' PSD Power')
        fig2.suptitle(title+' Piezo Power')
        pp.savefig(fig1)
        pp.savefig(fig2)
        del(LBD)
    pp.close()
