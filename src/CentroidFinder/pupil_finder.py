import Graffity
import matplotlib.pyplot as pyplot
import numpy
import scipy
import glob
import PIL
import images2gif
import io

fig = pyplot.figure(0)
fig.clear()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

# Change this directory to the location of the .TIF files on your system
directory= '/home/deen/Data/GRAVITY/Derotator/derotator/'
# Change this to match the structure of the files
files = glob.glob(directory+'pupil*001.TIF')

# These are the (x,y) coordinates and radius of the pupil cutout
xinit = 156
yinit = 87
width = 60
# Zoom Factor - how finely you want to sample the grid
ZF = 3.0

# Want to check how the Pupil Image looks like?  Set to True
checkPupil = False

image = Graffity.FLIRCamImage(files[0])
pupilImage = image.zoomIn(image.extractCutout(xinit, yinit, width, chopTop=True), ZF)

if checkPupil:
    ax.matshow(pupilImage)
    fig.show()
    raw_input()
    ax.clear()

# Measures the center of mass of the Pupil Image
center = scipy.ndimage.measurements.center_of_mass(pupilImage)
center_x = center[0]/ZF - width
center_y = center[1]/ZF - width


# GUESSES for the rotation axis (in pixels)
XGUESS = 167
YGUESS = 87

x = []
y = []
sigma = []
angle = []
cutouts = []


# Start the loop
for df in files:
    image = Graffity.FLIRCamImage(df)        # Read in the next file
    
    # parse the angle from the filename
    angle.append(float(df.split('\\')[-1].split()[1].split('deg')[0]))
    
    # Finds the center of the image by cross-correlation routine
    centroid = image.findPupilCenter(XGUESS, YGUESS, zoomFactor = ZF, pupilImage=pupilImage)
    x.append(centroid[0]+center_x)
    y.append(centroid[1]+center_y)
    print x[-1], y[-1]
    cutouts.append(image.imdata)

ax.clear()
x = numpy.array(x)
y = numpy.array(y)
angle = numpy.array(angle)
order = numpy.argsort(angle)
ax.plot(x[order], y[order])


frames = []
buf = []
outfile = open('Pupil_Runout.txt', 'w')

for i in order:
    ax.clear()
    ax.matshow(cutouts[i])
    ax.plot(x[order], y[order], color = 'y', lw=4.0)
    ax.scatter(x[i], y[i], color = 'k', s=85.0)
    ax.set_xbound(lower=95, upper=230)
    ax.set_ybound(lower=20, upper=160)
    ax.text(100.0, 25.0, 'Angle = %d' % angle[i], fontsize=16, color = 'y')
    ax.text(100.0, 30.0, 'X = %.2f' % x[i], fontsize=16, color = 'y')
    ax.text(100.0, 35.0, 'Y = %.2f' % y[i], fontsize=16, color = 'y')
    buf.append(io.BytesIO())
    fig.savefig(buf[-1], format='png')
    buf[-1].seek(0)
    print angle[i], x[i], y[i]
    outfile.write("%d %.2f %.2f\n" % (angle[i], x[i], y[i]))
    frames.append(PIL.Image.open(buf[-1]))

images2gif.writeGif('Pupil_Runout.gif', frames, duration=0.5)


