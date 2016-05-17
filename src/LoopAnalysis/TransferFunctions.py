import Graffity
import matplotlib.pyplot as pyplot
import glob
import numpy
import scipy.fftpack as fftpack

def G(s):
    T = 0.004*6.28
    tau = 6.8e-4*6.28
    retval = 0.5*numpy.exp(-tau*s)*(1.0-numpy.exp(-T*s))/(T**2.0*s**2.0)
    return retval

pyplot.rcParams["font.size"] = 32.0

datadir = '/home/deen/Data/GRAVITY/MATLAB_DATA/2016-03-01/IMTRACKING-165152/'
LoopData = glob.glob(datadir+'CIAO_LOOP_*.fits')
LoopData.sort()

S2M = glob.glob(datadir+'RecnOptimiser.S2M*.fits')
S2M.sort()

LD = []

for data, s2m in zip(LoopData, S2M):
    LD.append(Graffity.CircularBuffer(data, s2m))
    LD[-1].calculatePSD()


fig = pyplot.figure(0)
fig.clear()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

old = LD[0]
freq = LD[0].freq

pos = freq > 0

#S = fftpack.fft(old.time)
laplace_G = G(freq[pos])

#theoretical = laplace_G/(laplace_G/(1.0+laplace_G)*numpy.exp(-0.002*freq[pos]/2.0)*
#        (1.0+laplace_G))

H = laplace_G/(1.0+laplace_G)
WFS = laplace_G/((1.0+laplace_G)*numpy.exp(-0.004*6.28*freq[pos]/2.0))

theoretical = laplace_G/(WFS*(1.0+laplace_G))



for cld, old in zip(LD[-1].FFT[0:4], LD[0].FFT[0:4]):
    ax.plot(freq, cld/old, lw=2.0)

ax.plot(freq[pos], theoretical, color = 'k', lw = 3.0 )
#ax.plot(freq[pos], laplace_G, color = 'g', lw = 3.0 )
#ax.plot(freq[pos], WFS, color = 'b', lw = 3.0 )
#ax.plot(freq[pos], H, color = 'r', lw = 3.0 )
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xbound(0.8, 150.0)
ax.set_ybound(1e-4, 1e2)
ax.set_title("Modes 1-4")
fig.show()
fig.savefig("TF_CD.png")
