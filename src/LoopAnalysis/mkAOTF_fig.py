from matplotlib import pyplot
import scipy
import numpy

fig = pyplot.figure(0)
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
ax.clear()

df = open('/home/cdeen/Data/CIAO/AOTF_output.txt', 'r')

freq = []
modes = []
theoretical = []
for i in range(60):
    modes.append([])

for line in df.readlines():
    l = line.split()
    freq.append(float(l[0]))
    for i in range(60):
        modes[i].append(float(l[i+1]))
    theoretical.append(float(l[-1]))

freq = numpy.array(freq)
theoretical = numpy.array(theoretical)
for i in range(60):
    modes[i] = numpy.array(modes[i])


for i in range(45):
    #ax.scatter(freq, modes[i], lw=0.5)
    ax.plot(freq, modes[i], lw=0.25)

ax.plot(freq, theoretical, lw=2.0, color = 'k')
ax.loglog()
ax.set_xbound(1.0, 499.9678/2.0)
#ax.grid(True, which='both')
ax.set_title("On-Sky Transfer Functions, gain=0.5")
ax.set_xlabel("Frequency (Hz)")
ax.set_ylabel("Rejection")
fig.show()
fig.savefig("AOTF.png")

