import scipy
import numpy
import matplotlib.pyplot as pyplot
from matplotlib.collections import PatchCollection

class Derotator( object ):
    def __init__(self, misalignment):
        self.Z = 10.0   # mm
        self.Y = 17.3   # mm
        self.origin = Point(0.0, 0.0, 0.0)
        self.theta = numpy.arctan2(self.Z, self.Y)+misalignment
        self.beta = numpy.pi/4.0 + self.theta/2.0
        #print numpy.rad2deg(self.beta)

        self.m1Point = Point(0.0, 0.0, self.Z)
        self.m2Point = Point(0.0, self.Y, 0.0)
        self.m3Point = Point(0.0, 0.0, -self.Z)

        self.m1Normal = numpy.array([0.0, numpy.sin(self.beta), numpy.cos(self.beta)])
        self.m2Normal = numpy.array([0.0, -1.0, 0.0])
        self.m3Normal = numpy.array([0.0, numpy.sin(self.beta), -numpy.cos(self.beta)])

        self.makeMirrors()

    def makeMirrors(self):
        self.mirror1 = Plane(self.m1Point, self.m1Normal)
        self.mirror2 = Plane(self.m2Point, self.m2Normal)
        self.mirror3 = Plane(self.m3Point, self.m3Normal)

    def rotate(self, angle):
        self.mirror1.rotate(self.origin, angle, numpy.array([0.0, 0.0, 1.0]))
        self.mirror2.rotate(self.origin, angle, numpy.array([0.0, 0.0, 1.0]))
        self.mirror3.rotate(self.origin, angle, numpy.array([0.0, 0.0, 1.0]))

    def propogate(self, line):
        r1 = self.mirror1.reflection(line)
        r2 = self.mirror2.reflection(r1)
        r3 = self.mirror3.reflection(r2)
        return r3

    def translate(self, dx, dy, dz):
        shift = Point(dx, dy, dz)
        self.origin = self.origin+shift
        self.m1Point = self.m1Point+shift
        self.m2Point = self.m2Point+shift
        self.m3Point = self.m3Point+shift
        self.makeMirrors()

    def tiptilt(self, angle, axis):
        ux = axis[0]
        uy = axis[1]
        uz = axis[2]
        rotation_matrix = numpy.array([
            [numpy.cos(angle) + ux**2*(1.0-numpy.cos(angle)), 
                ux*uy*(1.0-numpy.cos(angle)) - uz*numpy.sin(angle),
                ux*uz*(1.0-numpy.cos(angle)) + uy*numpy.sin(angle)],
            [ux*uy*(1.0-numpy.cos(angle))+uz*numpy.sin(angle),
                numpy.cos(angle)+uy**2*(1.0-numpy.cos(angle)),
                uy*uz*(1.0-numpy.cos(angle))-ux*numpy.sin(angle)],
            [uz*ux*(1.0-numpy.cos(angle))-uy*numpy.sin(angle),
                uz*uy*(1.0-numpy.cos(angle))+ux*numpy.sin(angle),
                numpy.cos(angle)+uz**2*(1.0-numpy.cos(angle))]])
        v1 = self.origin- self.m1Point
        v1 = numpy.array([v1.x, v1.y, v1.z])
        new_vector = numpy.dot(rotation_matrix, v1)
        self.m1Normal = numpy.dot(rotation_matrix, self.m1Normal)
        v1 = Point(new_vector[0], new_vector[1], new_vector[2])
        self.m1Point = self.origin + v1
        v2 = self.origin- self.m2Point
        v2 = numpy.array([v2.x, v2.y, v2.z])
        new_vector = numpy.dot(rotation_matrix, v2)
        self.m2Normal = numpy.dot(rotation_matrix, self.m2Normal)
        v2 = Point(new_vector[0], new_vector[1], new_vector[2])
        self.m2Point = self.origin + v2
        v3 = self.origin- self.m3Point
        v3 = numpy.array([v3.x, v3.y, v3.z])
        new_vector = numpy.dot(rotation_matrix, v3)
        self.m3Normal = numpy.dot(rotation_matrix, self.m3Normal)
        v3 = Point(new_vector[0], new_vector[1], new_vector[2])
        self.m3Point = self.origin + v3

        self.makeMirrors()


class Point( object ):
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def __add__(self, other):
        return Point(self.x+other.x, self.y+other.y, self.z+other.z)

    def __sub__(self, other):
        return Point(self.x-other.x, self.y-other.y, self.z-other.z)

    def __repr__(self):
        return "x: "+str(self.x)+" y: "+str(self.y)+" z: "+str(self.z)
    def __str__(self):
        return "x: "+str(self.x)+" y: "+str(self.y)+" z: "+str(self.z)

class Line( object):
    def __init__(self, p, slope):
        norm = numpy.sqrt(slope[0]**2.0 + slope[1]**2.0 + slope[2]**2.0)
        self.slope = slope/norm
        self.a = slope[0]/norm
        self.b = slope[1]/norm
        self.c = slope[2]/norm
        self.p = p

    def traverse(self, t):
        newx = self.p.x + self.a*t
        newy = self.p.y + self.b*t
        newz = self.p.z + self.c*t
        newpoint = Point(newx, newy, newz)
        return newpoint

class Plane( object ):
    def __init__(self, p, normal):
        self.p = p
        self.normal = normal/numpy.sqrt(numpy.sum(numpy.square(normal)))
        self.calculatePlaneEqn()

    def rotate(self, p, angle, axis):
        vector = numpy.array([p.x-self.p.x, p.y-self.p.y, p.z-self.p.z])
        ux = axis[0]
        uy = axis[1]
        uz = axis[2]
        rotation_matrix = numpy.array([
            [numpy.cos(angle) + ux**2*(1.0-numpy.cos(angle)), 
                ux*uy*(1.0-numpy.cos(angle)) - uz*numpy.sin(angle),
                ux*uz*(1.0-numpy.cos(angle)) + uy*numpy.sin(angle)],
            [ux*uy*(1.0-numpy.cos(angle))+uz*numpy.sin(angle),
                numpy.cos(angle)+uy**2*(1.0-numpy.cos(angle)),
                uy*uz*(1.0-numpy.cos(angle))-ux*numpy.sin(angle)],
            [uz*ux*(1.0-numpy.cos(angle))-uy*numpy.sin(angle),
                uz*uy*(1.0-numpy.cos(angle))+ux*numpy.sin(angle),
                numpy.cos(angle)+uz**2*(1.0-numpy.cos(angle))]])
        new_vector = numpy.dot(rotation_matrix, vector)
        self.p = Point(p.x-new_vector[0], p.y-new_vector[1], p.z-new_vector[2])
        self.normal = numpy.dot(rotation_matrix, self.normal)
        self.calculatePlaneEqn()

    def calculatePlaneEqn(self):
        Ax = self.normal[0]
        Ay = self.normal[1]
        Az = self.normal[2]
        if (numpy.abs(Ax) < 1e-5) & (numpy.abs(Ay) < 1e-5):
            self.coeffs = numpy.array([0.0, 0.0, 1.0/self.p.z])
        elif (numpy.abs(Ax) < 1e-5) & (numpy.abs(Az) < 1e-5 ):
            self.coeffs = numpy.array([0.0, 1.0/self.p.y, 0.0])
        elif (numpy.abs(Ay) < 1e-5) & (numpy.abs(Az) < 1e-5):
            self.coeffs = numpy.array([1.0/self.p.x, 0.0, 0.0])
        elif (numpy.abs(Ax) < 1e-5):
            p2 = Point(self.p.x, self.p.y+1.0, self.p.z-Ay/Az)
            X = numpy.array([
                [self.p.y, self.p.z],
                [p2.y, p2.z]])
            Y = numpy.array([1.0, 1.0])
            BC = numpy.linalg.solve(X, Y)
            self.coeffs = numpy.array([0.0, BC[0], BC[1]])
        elif (numpy.abs(Ay) < 1e-5):
            p2 = Point(self.p.x+1.0, self.p.y, self.p.z-Ax/Az)
            X = numpy.array([
                [self.p.x, self.p.z],
                [p2.x, p2.z]])
            Y = numpy.array([1.0, 1.0])
            AC = numpy.linalg.solve(X, Y)
            self.coeffs = numpy.array([AC[0], 0.0, AC[1]])
        elif (numpy.abs(Az) < 1e-5):
            p2 = Point(self.p.x-Ay/Ax, self.p.y+1.0, self.p.z)
            X = numpy.array([
                [self.p.x, self.p.y],
                [p2.x, p2.y]])
            Y = numpy.array([1.0, 1.0])
            AB = numpy.linalg.solve(X,Y)
            self.coeffs = numpy.array([AB[0], AB[1], 0.0])
        else:
            p2 = Point(self.p.x, self.p.y+1.0, self.p.z-Ay/Az)
            p3 = Point(self.p.x+1.0, self.p.y, self.p.z-Ax/Az)
            p4 = Point(self.p.x-Ay/Ax, self.p.y+1.0, self.p.z)
            X = numpy.array([
                [p2.x, p2.y, p2.z],
                [p3.x, p3.y, p3.z],
                [p4.x, p4.y, p4.z]])
            Y = numpy.array([1.0, 1.0, 1.0])
            self.coeffs = numpy.linalg.solve(X, Y)

    def intersection(self, line):
        dotproduct = numpy.dot(self.normal, line.slope)
        if dotproduct == 0.0:
            return False
        else:
            D = (self.coeffs[0]*line.p.x + self.coeffs[1]*line.p.y +
                    self.coeffs[2]*line.p.z)
            t = (1.0-D)/(self.coeffs[0]*line.a + self.coeffs[1]*line.b +
                    self.coeffs[2]*line.c)
            return line.traverse(t)

    def reflection(self, line):
        reflection = self.intersection(line)
        if reflection:
            dot = numpy.dot(self.normal, line.slope)
            newx = line.a - 2*dot*self.normal[0]
            newy = line.b - 2*dot*self.normal[1]
            newz = line.c - 2*dot*self.normal[2]
            newLine = Line(reflection, [newx, newy, newz])
            return newLine
        else:
            print "Error! Line does not intersect plane!"
            return False

        
focalLength = 1370.0   #mm

pathLength = 20.57     #mm introduced by Derotator
fudgeFactor = 0.6
PM = 1340.0

pupilDistance = 14.7

origin = Point(0.0, 0.0, 0.0)
p1 = Point(0.0, 0.0, PM)
p2 = Point(0.00, 50.0, PM)
p3 = Point(50.0, 0.0, PM)
p4 = Point(0.00, -50.0, PM)
p5 = Point(-50.0, 0.0, PM)

derot = Derotator(misalignment=0.000)

beam_angle = numpy.arctan(50.0/1370.0)

l0 = Line(origin, [0.0, 0.0, -1.0])
l1 = Line(p1, [0.0, 0.0, -1.0])
l2 = Line(p2, [0.0, -beam_angle, -1.0])
l3 = Line(p3, [-beam_angle, 0.0, -1.0])
l4 = Line(p4, [0.0, beam_angle, -1.0])
l5 = Line(p5, [beam_angle, 0.0, -1.0])


focalPlane = Plane(Point(0.0, 0.0, PM - focalLength-fudgeFactor), numpy.array([0.0, 0.0, 1.0]))
focalPlane_DR = Plane(Point(0.0, 0.0, PM-focalLength-pathLength), numpy.array([0.0, 0.0, 1.0]))

pupilPlane = Plane(Point(0.0, 0.0, PM-focalLength-fudgeFactor+pupilDistance), numpy.array([0.0, 0.0, 1.0]))
pupilPlane_DR = Plane(Point(0.0, 0.0, PM-focalLength-pathLength+pupilDistance), numpy.array([0.0, 0.0, 1.0]))

nsteps = 10
dangle = 2.0*numpy.pi/nsteps

angle = []
a = []
b = []
c = []
d = []
e = []

theta = 0.0

derot.tiptilt(0.00, numpy.array([1.0, 0.0, 0.0]))
derot.translate(0.00, 0.0, 0.00)

#derot.rotate(numpy.deg2rad(45.0))
#spot = detectorPlane.intersection(derot.propogate(l2))
#print asdf

for i in range(nsteps)[0:5]:
    angle.append(theta)
    #spot = detectorPlane.intersection(derot.propogate(l1))
    spot = focalPlane_DR.intersection(derot.propogate(l1))
    a.append([spot.x, spot.y])
    #spot = detectorPlane.intersection(derot.propogate(l2))
    spot = pupilPlane_DR.intersection(derot.propogate(l2))
    b.append([spot.x, spot.y])
    #spot = detectorPlane.intersection(derot.propogate(l3))
    spot = pupilPlane_DR.intersection(derot.propogate(l3))
    c.append([spot.x, spot.y])
    #spot = detectorPlane.intersection(derot.propogate(l4))
    spot = pupilPlane_DR.intersection(derot.propogate(l4))
    d.append([spot.x, spot.y])
    #spot = detectorPlane.intersection(derot.propogate(l5))
    spot = pupilPlane_DR.intersection(derot.propogate(l5))
    e.append([spot.x, spot.y])

    derot.rotate(dangle)
    theta += dangle

    

a = numpy.array(a)
b = numpy.array(b)
c = numpy.array(c)
d = numpy.array(d)
e = numpy.array(e)

pyplot.ion()
fig = pyplot.figure(0)
fig.clear()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

"""
s1 = focalPlane.intersection(l1)
s2 = focalPlane.intersection(l2)
s3 = focalPlane.intersection(l3)
s4 = focalPlane.intersection(l4)
s5 = focalPlane.intersection(l5)
ax.scatter([s1.x, s2.x, s3.x, s4.x, s5.x],[s1.y, s2.y, s3.y, s4.y, s5.y], color='r')

s1 = focalPlane_DR.intersection(derot.propogate(l1))
s2 = focalPlane_DR.intersection(derot.propogate(l2))
s3 = focalPlane_DR.intersection(derot.propogate(l3))
s4 = focalPlane_DR.intersection(derot.propogate(l4))
s5 = focalPlane_DR.intersection(derot.propogate(l5))
ax.scatter([s1.x, s2.x, s3.x, s4.x, s5.x],[s1.y, s2.y, s3.y, s4.y, s5.y], color='b')

s1 = pupilPlane.intersection(l1)
s2 = pupilPlane.intersection(l2)
s3 = pupilPlane.intersection(l3)
s4 = pupilPlane.intersection(l4)
s5 = pupilPlane.intersection(l5)
ax.scatter([s1.x, s2.x, s3.x, s4.x, s5.x],[s1.y, s2.y, s3.y, s4.y, s5.y], color='g')

s1 = pupilPlane_DR.intersection(derot.propogate(l1))
s2 = pupilPlane_DR.intersection(derot.propogate(l2))
s3 = pupilPlane_DR.intersection(derot.propogate(l3))
s4 = pupilPlane_DR.intersection(derot.propogate(l4))
s5 = pupilPlane_DR.intersection(derot.propogate(l5))
ax.scatter([s1.x, s2.x, s3.x, s4.x, s5.x],[s1.y, s2.y, s3.y, s4.y, s5.y], color='m')

fig.show()

print asdf
"""
for spots in zip(a, b, c, d, e):
    sp = numpy.array(spots[1:])
    center = numpy.mean(sp, axis=0)
    junk = pyplot.Circle(center, 1.2, lw=10.0, ec='r', fc='b', alpha=0.5, fill=False)
    ax.add_artist(junk)
    ax.scatter(center[0], center[1], color='r')
    ax.scatter(spots[0][0], spots[0][1], color='b')
    print spots[0] - center
    #ax.plot(sp[:,0], sp[:,1], marker='o')


ax.set_xlim(-1.5, 1.5)
ax.set_ylim(-1.5, 1.5)
fig.show()
"""
ax.plot(a[:,0], a[:,1], c='r', marker='o')
ax.plot(b[:,0], b[:,1], c='b', marker='x')
ax.plot(c[:,0], c[:,1], c='g', marker='+')
ax.plot(d[:,0], d[:,1], c='k', marker='.')
#ax.plot(e[:,0], e[:,1], c='m', marker='x')
"""
