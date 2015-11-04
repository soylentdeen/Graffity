import Graffity
import numpy
import scipy
import matplotlib.pyplot as pyplot

wave = 632.8

ciao = Graffity.WFS(wavelength=1800.0)

var = numpy.array([False, False, True, True, True])
offsets = []
x = 0
for v in var:
    if v:
        offsets.append(x)
        x+= 1
    else:
        offsets.append(0)

zern = [0.0, 0.0, 0.0, 0.0, 0.0]
pupil = [0.0, 0.0]
actPoke = numpy.zeros(60, dtype=numpy.float32)
derotAngle = 0.00
clockingAngle = 0.0

# Take the flat-wavefront image
ciao.setupInstrument(zern, pupil, actPoke, derotAngle, clockingAngle)
ciao.expose()

if var[0]:
    f0 = pyplot.figure(0)
    f0.clear()
    ax0 = f0.add_axes([0.1, 0.1, 0.8, 0.8])
if var[1]:
    f1 = pyplot.figure(1)
    f1.clear()
    ax1 = f1.add_axes([0.1, 0.1, 0.8, 0.8])
if var[2]:
    f2 = pyplot.figure(2)
    f2.clear()
    ax2 = f2.add_axes([0.1, 0.1, 0.8, 0.8])
if var[3]:
    f3 = pyplot.figure(3)
    f3.clear()
    ax3 = f3.add_axes([0.1, 0.1, 0.8, 0.8])
if var[4]:
    f4 = pyplot.figure(4)
    f4.clear()
    ax4 = f4.add_axes([0.1, 0.1, 0.8, 0.8])
f5 = pyplot.figure(5)
f5.clear()
f6 = pyplot.figure(6)
f6.clear()

ax5 = f5.add_axes([0.1, 0.1, 0.8, 0.8])
ax6 = f6.add_axes([0.1, 0.1, 0.8, 0.8])

wferror = numpy.linspace(-2.0*wave, 2.0*wave, num=14)

clockingAngle = 0.00
for rms in wferror:
    print rms
    if var[0]:
        zern = [rms, 0.0, 0.0, 0.0, 0.0]
        ciao.setupInstrument(zern, pupil, actPoke, derotAngle, clockingAngle)
        ciao.expose()
    if var[1]:
        zern = [0.0, rms, 0.0, 0.0, 0.0]
        ciao.setupInstrument(zern, pupil, actPoke, derotAngle, clockingAngle)
        ciao.expose()
    if var[2]:
        zern = [0.0, 0.0, rms, 0.0, 0.0]
        ciao.setupInstrument(zern, pupil, actPoke, derotAngle, clockingAngle)
        ciao.expose()
    if var[3]:
        zern = [0.0, 0.0, 0.0, rms, 0.0]
        ciao.setupInstrument(zern, pupil, actPoke, derotAngle, clockingAngle)
        ciao.expose()
    if var[4]:
        zern = [0.0, 0.0, 0.0, 0.0, rms]
        ciao.setupInstrument(zern, pupil, actPoke, derotAngle, clockingAngle)
        ciao.expose()

centroids = numpy.array(ciao.centroids)
nvars = len(var[var==True])
flat = centroids[0]
if var[0]:
    tip = centroids[[i*nvars+offsets[0]+1 for i in range(len(wferror))]]-flat
if var[1]:
    tilt = centroids[[i*nvars+offsets[1]+1 for i in range(len(wferror))]]-flat
if var[2]:
    focus = centroids[[i*nvars+offsets[2]+1 for i in range(len(wferror))]]-flat
if var[3]:
    astig1 = centroids[[i*nvars+offsets[3]+1 for i in range(len(wferror))]]-flat
if var[4]:
    astig2 = centroids[[i*nvars+offsets[4]+1 for i in range(len(wferror))]]-flat


colorMap = pyplot.get_cmap()
colors = [colorMap(i) for i in numpy.linspace(0, 1, len(wferror))]

subapnum = range(68)

rms_x = []
rms_y = []
max_x = []
max_y = []

for i in range(len(wferror)):
    rx = []
    ry = []
    mx = []
    my = []
    if var[0]:
        ax0.plot(subapnum, tip[:,:,0][i], color=colors[i], marker='o')
        ax0.plot(subapnum, tip[:,:,1][i], color=colors[i], marker='+')
        rx.append(numpy.sqrt(numpy.average(tip[:,:,0][i]**2.0)))
        ry.append(numpy.sqrt(numpy.average(tip[:,:,1][i]**2.0)))
        mx.append(numpy.max(numpy.abs(tip[:,:,0][i])))
        my.append(numpy.max(numpy.abs(tip[:,:,1][i])))
    if var[1]:
        ax1.plot(subapnum, tilt[:,:,0][i], color=colors[i], marker='o')
        ax1.plot(subapnum, tilt[:,:,1][i], color=colors[i], marker='+')
        rx.append(numpy.sqrt(numpy.average(tilt[:,:,0][i]**2.0)))
        ry.append(numpy.sqrt(numpy.average(tilt[:,:,1][i]**2.0)))
        mx.append(numpy.max(numpy.abs(tilt[:,:,0][i])))
        my.append(numpy.max(numpy.abs(tilt[:,:,1][i])))
    if var[2]:
        ax2.plot(subapnum, focus[:,:,0][i], color=colors[i], marker='o')
        ax2.plot(subapnum, focus[:,:,1][i], color=colors[i], marker='+')
        rx.append(numpy.sqrt(numpy.average(focus[:,:,0][i]**2.0)))
        ry.append(numpy.sqrt(numpy.average(focus[:,:,1][i]**2.0)))
        mx.append(numpy.max(numpy.abs(focus[:,:,0][i])))
        my.append(numpy.max(numpy.abs(focus[:,:,1][i])))
    if var[3]:
        ax3.plot(subapnum, astig1[:,:,0][i], color=colors[i], marker='o')
        ax3.plot(subapnum, astig1[:,:,1][i], color=colors[i], marker='+')
        rx.append(numpy.sqrt(numpy.average(astig1[:,:,0][i]**2.0)))
        ry.append(numpy.sqrt(numpy.average(astig1[:,:,1][i]**2.0)))
        mx.append(numpy.max(numpy.abs(astig1[:,:,0][i])))
        my.append(numpy.max(numpy.abs(astig1[:,:,1][i])))
    if var[4]:
        ax4.plot(subapnum, astig2[:,:,0][i], color=colors[i], marker='o')
        ax4.plot(subapnum, astig2[:,:,1][i], color=colors[i], marker='+')
        rx.append(numpy.sqrt(numpy.average(astig2[:,:,0][i]**2.0)))
        ry.append(numpy.sqrt(numpy.average(astig2[:,:,1][i]**2.0)))
        mx.append(numpy.max(numpy.abs(astig2[:,:,0][i])))
        my.append(numpy.max(numpy.abs(astig2[:,:,1][i])))
    rms_x.append(rx)
    rms_y.append(ry)
    max_x.append(mx)
    max_y.append(my)


rms_x = numpy.array(rms_x).transpose()
rms_y = numpy.array(rms_y).transpose()
max_x = numpy.array(max_x).transpose()
max_y = numpy.array(max_y).transpose()

labels = []
lines = []
if var[0]:
    lines.append(ax5.plot(wferror, max_x[offsets[0]], color = 'b', marker = 'o')[0])
    ax6.plot(wferror, max_y[offsets[0]], color = 'b', marker = 'o')
    labels.append["Tip"]
    #ax5.plot(wferror, rms_x[0], color = 'b', marker = '+')
if var[1]:
    lines.append(ax5.plot(wferror, max_x[offsets[1]], color = 'g', marker = 'o')[0])
    ax6.plot(wferror, max_y[offsets[1]], color = 'g', marker = 'o')
    labels.append["Tilt"]
    #ax5.plot(wferror, rms_x[1], color = 'g', marker = '+')
if var[2]:
    lines.append(ax5.plot(wferror, max_x[offsets[2]], color = 'r', marker = 'o')[0])
    ax6.plot(wferror, max_y[offsets[2]], color = 'r', marker = 'o')
    labels.append("Focus")
    #ax5.plot(wferror, rms_x[2], color = 'r', marker = '+')
if var[3]:
    lines.append(ax5.plot(wferror, max_x[offsets[3]], color = 'c', marker = 'o')[0])
    ax6.plot(wferror, max_y[offsets[3]], color = 'c', marker = 'o')
    labels.append("Astig1")
    #ax5.plot(wferror, rms_x[3], color = 'c', marker = '+')
if var[4]:
    lines.append(ax5.plot(wferror, max_x[offsets[4]], color = 'm', marker = 'o')[0])
    ax6.plot(wferror, max_y[offsets[4]], color = 'm', marker = 'o')
    labels.append("Astig2")
    #ax5.plot(wferror, rms_x[4], color = 'm', marker = '+')

ax5.set_xlabel("RMS Wavefront Error (nm)")
ax5.set_ylabel("Maximum X Slope (pixels)")

f5.legend(lines, labels)
ax5.set_title('X Slopes')

ax6.set_xlabel("RMS Wavefront Error (nm)")
ax6.set_ylabel("Maximum Y Slope (pixels)")

f6.legend(lines, labels)
ax6.set_title('Y Slopes')

if var[0]:
    ax0.set_xlabel("Subaperture Number")
    ax0.set_ylabel("Slopes (Pixels)")
    ax0.set_title('Tip')
    f0.show()
    f0.savefig('tip.png')
if var[1]:
    ax1.set_xlabel("Subaperture Number")
    ax1.set_ylabel("Slopes (Pixels)")
    ax1.set_title('Tilt')
    f1.show()
    f1.savefig('tilt.png')
if var[2]:
    ax2.set_xlabel("Subaperture Number")
    ax2.set_ylabel("Slopes (Pixels)")
    ax2.set_title('Focus')
    f2.show()
    f2.savefig('focus.png')
if var[3]:
    ax3.set_xlabel("Subaperture Number")
    ax3.set_ylabel("Slopes (Pixels)")
    ax3.set_title('Oblique Astigmatism')
    f3.show()
    f3.savefig('ObliqAstig.png')
if var[4]:
    ax4.set_xlabel("Subaperture Number")
    ax4.set_ylabel("Slopes (Pixels)")
    ax4.set_title('Vertical Astigmatism')
    f4.show()
    f4.savefig('VertAstig.png')

f5.show()
f5.savefig('Xerror.png')
f6.show()
f6.savefig('Yerror.png')
