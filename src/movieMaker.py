import scipy
import numpy
import matplotlib.pyplot as pyplot
import pyfits
import matplotlib.animation as manimation

FFMpegWriter = manimation.writers['ffmpeg']
metadata = dict(title='Move Test', artist ='Matplotlib', comment='Movie Support!')

writer = FFMpegWriter(fps=15, metadata=metadata)

data = pyfits.getdata('/home/deen/Data/GRAVITY/MATLAB_DATA/2015-10-02/CIAO_PIXELS_0001.fits')

frames = data.field(3)
fnum = data.field(0)

fig = pyplot.figure()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
l = ax.matshow(numpy.zeros((72, 72)),cmap = pyplot.cm.rainbow)
l.norm.vmin = 0.0
l.norm.vmax = 3000.0
#l.norm.vmin = -1.0
#l.norm.vmax = 5.0
ax.set_xlim(0, 72)
ax.set_ylim(0, 72)
ax.set_yticklabels([])
ax.set_xticklabels([])
ho = ax.text(75, 20, 'HO: OPEN')
tt = ax.text(75, 10, 'TT: OPEN')

def HO_open(i):
    if fnum[i] < 35840:
        return True
    if fnum[i] < 39687:
        return False
    if fnum[i] < 40886:
        return True
    return False

def TT_open(i):
    if fnum[i] < 36008:
        return True
    if fnum[i] < 39586:
        return False
    if fnum[i] < 41046:
        return True
    return False

with writer.saving(fig, "CIAO_LPE.mp4", 100):
    for i in range(len(frames)):
        if HO_open(i):
            ho.set_text('HO: OPEN')
        else:
            ho.set_text('HO: CLOSED')
        if TT_open(i):
            tt.set_text('TT: OPEN')
        else:
            tt.set_text('TT: CLOSED')
        #l.set_data(numpy.log10(frames[i].reshape((72, 72))+1.0))
        l.set_data(frames[i].reshape((72, 72)))
        writer.grab_frame()
