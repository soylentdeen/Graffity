#!/usr/bin/env python
# -*- coding: utf-8 -*-

#avg											110,429166666667	142,303939393939	39,4660202129861	17,303315566334							358,067727272727	145,667575757576	39,0507475838317	18,0599499387257							607,147575757576	146,828787878788	39,4196129733035	17,4218484639878							854,394393939394	146,564696969697	38,8185123201097	18,0525720203619


import pyfits as fits
import numpy as np
from matplotlib import pyplot as plt
from scipy import ndimage
from astropy.modeling import models, fitting, Fittable2DModel
import argparse

northern_stations=['G2', 'J3', 'J4', 'J5', 'J6']

approx_scale = {'AT' : 80., 'UT' : 18.}

roof_x = np.array([110.4, 358.1, 607.1, 854.4])
roof_y = np.array([142.3, 145.7, 146.8, 146.6])

roof_pos = np.array([38.49, 38.54, 38.76, 39.80])

def mean_rms(data, axis=None, mean_fcn=np.mean):
    m = mean_fcn(data, axis=axis)
    r = np.sqrt(mean_fcn(np.square(data-m), axis=axis))
    return m, r

def format_results(label, data, uncertainties, fmt_str):
    out = label
    for i in np.arange(len(data)):
        out += ' ' + fmt_str.format(data[i], uncertainties[i])
    return out

def fit_port(port, approx_PAp, roof_xp, roof_yp, rho, win, plot, array):
    
    approx_dx=rho*np.sin(approx_PAp*np.pi/180)/approx_scale[array];
    approx_dy=rho*np.cos(approx_PAp*np.pi/180)/approx_scale[array];
    
    xmax=int(roof_xp-0.5*approx_dx)
    ymax=int(roof_yp-0.5*approx_dy)

    
    thumb=port[ymax-win:ymax+win+1, xmax-win:xmax+win+1]

    y, x = np.mgrid[-win:win+1, -win:win+1]
    
    g_init=models.Gaussian2D(amplitude=thumb.max(), x_mean=0, y_mean=0, x_stddev=3., y_stddev=3., theta=0.)
    fit_g = fitting.LevMarLSQFitter()
    g = fit_g(g_init, x, y, thumb)
    
    if plot >= 2:
        plt.figure(figsize=(8, 2.5))
        plt.subplot(1, 3, 1)
        plt.imshow(thumb, origin='lower', interpolation='nearest', vmin=0, vmax=thumb.max())
        plt.title("Data")
        plt.subplot(1, 3, 2)
        plt.imshow(g(x, y), origin='lower', interpolation='nearest', vmin=0, vmax=thumb.max())
        plt.title("Model")
        plt.subplot(1, 3, 3)
        plt.imshow(thumb - g(x, y), origin='lower', interpolation='nearest', vmin=0, vmax=thumb.max())
        plt.title("Residual")
        plt.show()
    
    x1 = xmax + g.x_mean
    y1 = ymax + g.y_mean
    
        
    x2max=int(roof_xp+0.5*approx_dx)
    y2max=int(roof_yp+0.5*approx_dy)
    
    thumb2=port[y2max-win:y2max+win+1, x2max-win:x2max+win+1]
    
    g2_init=models.Gaussian2D(amplitude=thumb2.max(), x_mean=0, y_mean=0, x_stddev=g.x_stddev, y_stddev=g.y_stddev, theta=0.)
    g2_init.x_stddev.fixed=True
    g2_init.y_stddev.fixed=True
    g2 = fit_g(g2_init, x, y, thumb2)
    
    x2 = x2max+g2.x_mean
    y2 = y2max+g2.y_mean

    if plot >= 2:
        plt.figure(figsize=(8, 2.5))
        plt.subplot(1, 3, 1)
        plt.imshow(thumb2, origin='lower', interpolation='nearest', vmin=0, vmax=thumb2.max())
        plt.title("Data")
        plt.subplot(1, 3, 2)
        plt.imshow(g2(x, y), origin='lower', interpolation='nearest', vmin=0, vmax=thumb2.max())
        plt.title("Model")
        plt.subplot(1, 3, 3)
        plt.imshow(thumb2 - g2(x, y), origin='lower', interpolation='nearest', vmin=0, vmax=thumb2.max())
        plt.title("Residual")
        plt.show()
    
    return x1, y1, x2, y2

def do_all(data, dark, theta_in, rho_in, win, approx_PA, group, plot, array):
    
    ngroups=int(data.shape[0])/group
    if ngroups == 0:
        ngroups=1
    
    x1s=np.zeros((ngroups,4))
    y1s=np.zeros((ngroups,4))
    x2s=np.zeros((ngroups,4))
    y2s=np.zeros((ngroups,4))
    
    
    for g in np.arange(ngroups):
        avg_data=np.mean(data[g*group:(g+1)*group, :, :], axis=0)-dark
        for p in np.arange(4):
            port=avg_data[0:250, p*250:(p+1)*250]
            x1, y1, x2, y2 = fit_port(port, approx_PA[p], roof_x[p]-p*250, roof_y[p], rho_in, win, plot, array)
            x1s[g, p]=x1+p*250
            y1s[g, p]=y1
            x2s[g, p]=x2+p*250
            y2s[g, p]=y2
    
    return x1s, y1s, x2s, y2s

if (__name__ == "__main__"):
    parser = argparse.ArgumentParser(description='Measure star positions on acquisition camera.')
    parser.add_argument('fname', nargs='+')
    parser.add_argument("-w", "--window", help="full width of window for peak fitting", type=int)
    parser.add_argument("-g", "--group", help="number of frames to average together", type=int)
    parser.add_argument("-r", "--rho", help="separation of binary (mas)", type=float)
    parser.add_argument("-t", "--theta", help="position angle of binary on sky (East of North)", type=float)
    parser.add_argument("-d", "--dark_file", help="dark file", action='append')
    parser.add_argument("-D", "--dark_outfile", help="file to save preprocessed dark")
    parser.add_argument("-p", "--plot", help="what to plot. 0: nothing, 1: final fit results, 2: final results+individual fit residuals", type=int, default=0)
    parser.add_argument("-v", "--verbose", help="verbosity level.", type=int, default=0)

    args = parser.parse_args()
    args.window
    
    verbose=args.verbose

    if args.dark_file is None:
        dark=0.
    else:
        dark_file=args.dark_file
        infile=fits.open(dark_file.pop())
        prihdr=infile[0].header
        dark=infile['IMAGING_DATA_ACQ']
        darkhdr=dark.header
        dark=dark.data
        for f in dark_file:
            dark=np.append(dark, fits.open(f)['IMAGING_DATA_ACQ'].data, axis=0)
        NDIT=dark.shape[0]
        dark=np.median(dark, axis=0)
        if args.dark_outfile is not None:
            pri=fits.PrimaryHDU(header=prihdr)
            img=fits.ImageHDU(dark[None,:,:], header=darkhdr)
            pri.header.set('HIERARCH ESO DET1 NDIT', NDIT, 'number of DITs accumulated together')
            newfits=fits.HDUList([pri, img])
            newfits.writeto(args.dark_outfile)
        dark=dark[:250,:]
            
    
    x1s_mean=np.zeros((len(args.fname),4))
    x1s_rms=np.zeros((len(args.fname),4))
    x2s_mean=np.zeros((len(args.fname),4))
    x2s_rms=np.zeros((len(args.fname),4))
    y1s_mean=np.zeros((len(args.fname),4))
    y1s_rms=np.zeros((len(args.fname),4))
    y2s_mean=np.zeros((len(args.fname),4))
    y2s_rms=np.zeros((len(args.fname),4))
    for f in np.arange(len(args.fname)):
        fname=args.fname[f]
        hdulist=fits.open(fname)
        data=hdulist['IMAGING_DATA_ACQ'].data[:,:250,:]
        
        if args.group is None:
            group=data.shape[0]/4
        else:
            group=args.group
            
        if args.rho is None or args.theta is None:
            dx_in=hdulist[0].header["HIERARCH ESO INS SOBJ X"]
            dy_in=hdulist[0].header["HIERARCH ESO INS SOBJ Y"]
                
        if args.rho is None:
            rho_in=np.sqrt(dx_in*dx_in+dy_in*dy_in)
        else:
            rho_in=args.rho
                    
        enc=np.zeros(4)
        approx_PA=np.zeros(4)
        sta=dict()

        GVPORTS={7: 'GV1', 5: 'GV2', 3: 'GV3', 1:'GV4'}
        config=''
        array='AT'
        for t in ('4', '3', '2', '1'):
            port=GVPORTS[hdulist[0].header["HIERARCH ESO ISS CONF INPUT"+t]]
            tel=hdulist[0].header["HIERARCH ESO ISS CONF T"+t+"NAME"]
            if tel[0:2] == 'UT':
                array='UT'
            sta[port]=hdulist[0].header["HIERARCH ESO ISS CONF STATION"+t]
            config += ' {port}: {tel}/{sta}/{dl}'.format(port=port, tel=tel, sta=sta[port], dl=hdulist[0].header["HIERARCH ESO ISS CONF DL"+t])
        
        if args.theta is None:
            theta_in=np.arctan2(dx_in,dy_in)*180./np.pi
        else:
            theta_in=args.theta
                        
        if args.window is None:
            win=long(rho_in/2./approx_scale[array]/np.sqrt(2.))
        else:
            win=long(args.window)/2
        
#        if array == 'UT':
        if True:
            for p in np.arange(4):
                try:
                    rp=hdulist[0].header["HIERARCH ESO INS DROTOFF"+str(p+1)]
                except (KeyError):
                    rp=roof_pos[p]
                approx_PA[p] = 270.-rp
        else:
            for p in np.arange(4):
                enc[p]=hdulist[0].header["HIERARCH ESO INS DROT"+str(p+1)+" ENC"]
                approx_PA[p]=theta_in-enc[p]/200.;
                if sta['GV'+str(p+1)] in northern_stations:
                    approx_PA[p] += 180.

        if verbose:
            print("*******************************************************************************")
            print("File name     : {}".format(fname))
            print("FT target     : {}".format(hdulist[0].header["HIERARCH ESO FT ROBJ NAME"]))
            print("SC target     : {}".format(hdulist[0].header["HIERARCH ESO INS SOBJ NAME"]))
            print("Position angle: {}°".format(theta_in))
            print("Separation    : {} mas".format(rho_in))
            print("Configuration :{}".format(config))
  
        hdulist.close()

        x1s, y1s, x2s, y2s = do_all(data, dark, theta_in, rho_in, win, approx_PA, group, args.plot, array)
        
        
        x1s_mean[f,:], x1s_rms[f,:] = mean_rms(x1s, axis=0, mean_fcn=np.median)
        y1s_mean[f,:], y1s_rms[f,:] = mean_rms(y1s, axis=0, mean_fcn=np.median)
        x2s_mean[f,:], x2s_rms[f,:] = mean_rms(x2s, axis=0, mean_fcn=np.median)
        y2s_mean[f,:], y2s_rms[f,:] = mean_rms(y2s, axis=0, mean_fcn=np.median)
        
        if verbose:
            print(format_results('x1:', x1s_mean[f,:], x1s_rms[f,:], '{:8.3f}±{:.4f}'))
            print(format_results('y1:', y1s_mean[f,:], y1s_rms[f,:], '{:8.3f}±{:.4f}'))
            print(format_results('x2:', x2s_mean[f,:], x2s_rms[f,:], '{:8.3f}±{:.4f}'))
            print(format_results('y2:', y2s_mean[f,:], y2s_rms[f,:], '{:8.3f}±{:.4f}'))

        txt = '{} {:8s} {:8s} '.format(fname,
                                 hdulist[0].header["HIERARCH ESO FT ROBJ NAME"],
                                 hdulist[0].header["HIERARCH ESO INS SOBJ NAME"])
        for i in np.arange(4):
            txt += '{:8.3f} {:8.3f} {:8.3f} {:8.3f} '.format(x1s_mean[f,i], y1s_mean[f,i], x2s_mean[f,i], y2s_mean[f,i])
        print (txt)

