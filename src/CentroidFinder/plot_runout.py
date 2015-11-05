import scipy
import numpy
import matplotlib.pyplot as pyplot
import PIL
import images2gif
import io

fig = pyplot.figure(0)
fig.clear()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

pupilFile = open('Pupil_Runout.txt', 'r')
focusFile = open('Focus_Runout.txt', 'r')

pupilData = pupilFile.readlines()
focusData = focusFile.readlines()

pupilFile.close()
focusFile.close()

angle = []
pupil_x = []
pupil_y = []
focus_x = []
focus_y = []

for pupil, focus in zip(pupilData, focusData):
    pupdat = pupil.split()
    focdat = focus.split()
    angle.append(float(pupdat[0]))
    pupil_x.append(float(pupdat[1]))
    pupil_y.append(float(pupdat[2]))
    focus_x.append(float(focdat[1]))
    focus_y.append(float(focdat[2]))

angle = numpy.array(angle)
focus_x = numpy.array(focus_x)
focus_y = numpy.array(focus_y)
pupil_x = numpy.array(pupil_x)
pupil_y = numpy.array(pupil_y)

ax.plot(focus_x, focus_y)
ax.plot(pupil_x, pupil_y)

fig.savefig("pupil_and_focus.png")

ax.clear()

frames = []
buf = []

order = numpy.argsort(angle)
for i in order:
    ax.clear()
    ax.plot(focus_x, focus_y, color = 'b', lw=2.0)
    ax.plot(pupil_x, pupil_y, color = 'g', lw=2.0)
    ax.scatter(focus_x[i], focus_y[i], color = 'k', s=30.0)
    ax.scatter(pupil_x[i], pupil_y[i], color = 'k', s=30.0)
    ax.text(157, 82, 'Angle = %d' % angle[i])
    buf.append(io.BytesIO())
    fig.savefig(buf[-1], format='png')
    buf[-1].seek(0)
    frames.append(PIL.Image.open(buf[-1]))

images2gif.writeGif('FocusPupil_Runout.gif', frames, duration = 1.0)

