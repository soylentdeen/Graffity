import scipy
import numpy
from matplotlib import pyplot
import Graffity

fig = pyplot.figure(0)
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

PSF = Graffity.PSF()

PSF.generateOTF()

ax.matshow(PSF.getPSF())
fig.show()
