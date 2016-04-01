import scipy
import numpy
import matplotlib.pyplot as pyplot

pyplot.rc('font', size=19.0)
fig = pyplot.figure(0, figsize=(8,6))
fig.clear()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

lon = numpy.deg2rad(70.4)
lat = numpy.deg2rad(-24.6)

dec = numpy.deg2rad(-29.0)
#dec = numpy.deg2rad(-24.6)
#alpha = numpy.deg2rad(255.0)
alpha = numpy.deg2rad(0)

#h = t - ra
#sin(dec) = sin(alt)*sin(lat) + cos(alt)*cos(lat)*cos(az)
#cos(h) = (sin(alt) - sin(lat)*sin(dec))/(cos(dec)*cos(lat))


alt = []
az = []
RA = []
delta = []
#HA = numpy.linspace(-2.3, -1.5, num=1000)
HA = numpy.linspace(-numpy.pi, numpy.pi, num=1000)
delta_t = numpy.rad2deg(HA[-1]-HA[-2])/360.0*24.0*3600.0#/2.0
for t in HA:
    alt.append(numpy.arcsin(numpy.cos(dec)*numpy.cos(lat)*numpy.cos(t - alpha) +
        numpy.sin(lat)*numpy.sin(dec)))
    
    az.append(numpy.arccos((numpy.sin(dec) - numpy.sin(alt[-1])*numpy.sin(lat))/
        (numpy.cos(alt[-1])*numpy.cos(lat))))

    #RA.append(-numpy.rad2deg(alt[-1])+numpy.rad2deg(az[-1])+90)
    RA.append(numpy.rad2deg(az[-1])+90)
    
    if len(RA) == 1:
        delta.append(0.0)
    else:
        delta.append((RA[-1] - RA[-2])/delta_t)

ax.plot(numpy.rad2deg(HA)/360.0*24.0, RA)
#ax.plot(HA, RA)
#ax.plot(numpy.rad2deg(HA[1:])/360.0*24.0, numpy.array(delta[1:]), lw=2.0)
ax.set_xbound(lower=-12, upper=12)
#ax.plot(numpy.rad2deg(HA)/360.0*24.0, alt)
#ax.plot(numpy.rad2deg(HA)/360.0*24.0, az)
ax.set_ylabel("Field Rotation: Degrees per Second")
ax.set_xlabel("Hour Angle for the Galactic Center")
fig.show()
fig.savefig('galactic_center_field_rotation.png')
