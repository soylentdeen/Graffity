import scipy
import numpy
import matplotlib.pyplot as pyplot
from matplotlib.collections import PatchCollection
from scipy import optimize
from numpy import random

class Derotator( object ):
    def __init__(self, iMs = {}, tipAngle = 0.0,
            tipAxis=numpy.array([1.0, 0.0, 0.0]),
            tiltAngle = 0.0, tiltAxis=numpy.array([0.0, 1.0, 0.0]),
            x=0.0, y=0.0, z=0.0, ax=None):
        self.Z = 10.0   # mm
        self.Y = 17.3   # mm
        self.rotationAxis = numpy.array([0.0, 0.0, 1.0])
        self.origin = Point(0.0, 0.0, 0.0)
	self.theta = numpy.arctan2(self.Z, self.Y)
        self.beta = numpy.pi/4.0 + self.theta/2.0

        self.m1Point = Point(0.0, 0.0, self.Z+iMs["m1Point_Z"])
        self.m2Point = Point(0.0, self.Y+iMs["m2Point_Y"], 0.0)
        self.m3Point = Point(0.0, 0.0, -self.Z+iMs["m3Point_Z"])

        self.m1Normal = numpy.array([0.0, numpy.sin(self.beta), numpy.cos(self.beta)])
        self.m2Normal = numpy.array([0.0, -1.0, 0.0])
        self.m3Normal = numpy.array([0.0, numpy.sin(self.beta), -numpy.cos(self.beta)])

        self.translate(x, y, z)
        self.tiptilt(tipAngle, tipAxis, tiltAngle, tiltAxis)

        self.makeMirrors(iMs=iMs)

        self.generateLimacons(ax=ax)

    def makeMirrors(self, iMs={}):
        self.mirror1 = Plane(self.m1Point, self.m1Normal)
        self.mirror2 = Plane(self.m2Point, self.m2Normal)
        self.mirror3 = Plane(self.m3Point, self.m3Normal)

        self.mirror1.rotate(self.m1Point, iMs["m1Angle"], iMs["m1Axis"])
        self.mirror2.rotate(self.m2Point, iMs["m2Angle"], iMs["m2Axis"])
        self.mirror3.rotate(self.m3Point, iMs["m3Angle"], iMs["m3Axis"])

    def rotate(self, angle):
        self.mirror1.rotate(self.origin, angle, self.rotationAxis)
        self.mirror2.rotate(self.origin, angle, self.rotationAxis)
        self.mirror3.rotate(self.origin, angle, self.rotationAxis)

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

    def tiptilt(self, tipAngle, tipAxis, tiltAngle, tiltAxix):
        ux = tipAxis[0]
        uy = tipAxis[1]
        uz = tipAxis[2]
        tip_rotation_matrix = numpy.array([
            [numpy.cos(tipAngle) + ux**2*(1.0-numpy.cos(tipAngle)), 
                ux*uy*(1.0-numpy.cos(tipAngle)) - uz*numpy.sin(tipAngle),
                ux*uz*(1.0-numpy.cos(tipAngle)) + uy*numpy.sin(tipAngle)],
            [ux*uy*(1.0-numpy.cos(tipAngle))+uz*numpy.sin(tipAngle),
                numpy.cos(tipAngle)+uy**2*(1.0-numpy.cos(tipAngle)),
                uy*uz*(1.0-numpy.cos(tipAngle))-ux*numpy.sin(tipAngle)],
            [uz*ux*(1.0-numpy.cos(tipAngle))-uy*numpy.sin(tipAngle),
                uz*uy*(1.0-numpy.cos(tipAngle))+ux*numpy.sin(tipAngle),
                numpy.cos(tipAngle)+uz**2*(1.0-numpy.cos(tipAngle))]])
        ux = tiltAxis[0]
        uy = tiltAxis[1]
        uz = tiltAxis[2]
        tilt_rotation_matrix = numpy.array([
            [numpy.cos(tiltAngle) + ux**2*(1.0-numpy.cos(tiltAngle)), 
                ux*uy*(1.0-numpy.cos(tiltAngle)) - uz*numpy.sin(tiltAngle),
                ux*uz*(1.0-numpy.cos(tiltAngle)) + uy*numpy.sin(tiltAngle)],
            [ux*uy*(1.0-numpy.cos(tiltAngle))+uz*numpy.sin(tiltAngle),
                numpy.cos(tiltAngle)+uy**2*(1.0-numpy.cos(tiltAngle)),
                uy*uz*(1.0-numpy.cos(tiltAngle))-ux*numpy.sin(tiltAngle)],
            [uz*ux*(1.0-numpy.cos(tiltAngle))-uy*numpy.sin(tiltAngle),
                uz*uy*(1.0-numpy.cos(tiltAngle))+ux*numpy.sin(tiltAngle),
                numpy.cos(tiltAngle)+uz**2*(1.0-numpy.cos(tiltAngle))]])
        rotation_matrix = tip_rotation_matrix.dot(tilt_rotation_matrix)
        v1 = self.origin- self.m1Point
        v1 = numpy.array([v1.x, v1.y, v1.z])
        new_vector = numpy.dot(rotation_matrix, v1)
        self.m1Normal = numpy.dot(rotation_matrix, self.m1Normal)
        v1 = Point(new_vector[0], new_vector[1], new_vector[2])
        self.m1Point = self.origin - v1
        v2 = self.origin- self.m2Point
        v2 = numpy.array([v2.x, v2.y, v2.z])
        new_vector = numpy.dot(rotation_matrix, v2)
        self.m2Normal = numpy.dot(rotation_matrix, self.m2Normal)
        v2 = Point(new_vector[0], new_vector[1], new_vector[2])
        self.m2Point = self.origin - v2
        v3 = self.origin- self.m3Point
        v3 = numpy.array([v3.x, v3.y, v3.z])
        new_vector = numpy.dot(rotation_matrix, v3)
        self.m3Normal = numpy.dot(rotation_matrix, self.m3Normal)
        v3 = Point(new_vector[0], new_vector[1], new_vector[2])
        self.m3Point = self.origin - v3
        self.rotationAxis = numpy.dot(rotation_matrix, self.rotationAxis)

    def generateLimacons(self, ax=None):
        focalLength = 1370.0   #mm

        pathLength = 20.57     #mm introduced by Derotator
        fudgeFactor = 0.6
        PM = 1340.0

        pupilDistance = 14.7

        origin = Point(0.0, 0.0, 0.0)
        #p1 = Point(0.0, 0.0, PM)
        #p2 = Point(0.00, 50.0, PM)
        #p3 = Point(50.0, 0.0, PM)
        #p4 = Point(0.00, -50.0, PM)
        #p5 = Point(-50.0, 0.0, PM)
        l0 = Line(origin, [0.0, 0.0, -1.0])
        #lxp = Line(origin, [40.0, 0.0, focalLength])
        #lxm = Line(origin, [-40.0, 0.0, focalLength])
        #lyp = Line(origin, [0.0, 40.0, focalLength])
        #lym = Line(origin, [0.0, -40.0, focalLength])
        
        focalPlane = Plane(Point(0.0, 0.0, PM - focalLength),
            numpy.array([0.0, 0.0, 1.0]))

        pupilPlane = Plane(Point(0.0, 0.0, PM-focalLength+pupilDistance),
                numpy.array([0.0, 0.0, 1.0]))

        nsteps = 50
        dangle = 2.0*numpy.pi/nsteps

        angle = []
        focal_spots = []
        pupil_centroids = []
        th = []

        theta = 0.0

        for i in range(nsteps):
            angle.append(theta)
            focal_spot = focalPlane.intersection(self.propogate(l0))
            focal_spots.append([focal_spot.x, focal_spot.y])

            #pxp = pupilPlane.intersection(self.propogate(lxp))
            #pxm = pupilPlane.intersection(self.propogate(lxm))
            #pyp = pupilPlane.intersection(self.propogate(lyp))
            #pym = pupilPlane.intersection(self.propogate(lym))
            pupil_centroid = pupilPlane.intersection(self.propogate(l0))
            pupil_centroids.append([pupil_centroid.x, pupil_centroid.y])
            th.append(theta)

            self.rotate(dangle)
            theta += dangle
    
        self.focal_spots = numpy.array(focal_spots)
        self.pupil_centroids = numpy.array(pupil_centroids)
        self.theta = numpy.array(th)

        if ax!= None:
            ax.clear()
            ax.set_aspect('equal')
            ax.plot(self.focal_spots.T[0], self.focal_spots.T[1], color = 'r')
            ax.plot(self.pupil_centroids.T[0], self.pupil_centroids.T[1], color = 'b')
            ax.figure.show()
            raw_input()
        
    """
    def fitfunc(self, p, theta):
        anglex = theta + p[6]
	angley = theta + p[7]
        x = p[0] + numpy.abs(p[1])*numpy.cos(anglex) + \
                numpy.abs(p[2])/2.0*numpy.cos(2.0*anglex)
        y = p[3] + numpy.abs(p[4])*numpy.sin(angley) + \
                numpy.abs(p[5])/2.0*numpy.sin(2.0*angley)
	R = numpy.array([[numpy.cos(p[8]), -numpy.sin(p[8])],
		          [numpy.sin(p[8]), numpy.cos(p[8])]])
	I = numpy.array([[1.0, 0.0],[0.0, 1.0]])
	X = numpy.array([x, y])
	X0 = numpy.array([p[0], p[3]])

	rotated = R.dot(X)
	translated = (I - R).dot(X0)
        return rotated[0]+translated[0], rotated[1]+translated[1]
    #"""

    def fitfunc(self, p, theta):
        C = p[0]*(1.0 + 0.0j) + p[3]*(0.0 + 1.0j)
        A = p[1]*(1.0 + 0.0j) + p[4]*(0.0 + 1.0j)
        B = p[2]*(1.0 + 0.0j) + p[5]*(0.0 + 1.0j)
        #A = numpy.abs(p[1])*(1.0 + 0.0j) + numpy.abs(p[4])*(0.0 + 1.0j)
        #B = numpy.abs(p[2])*(1.0 + 0.0j) + numpy.abs(p[5])*(0.0 + 1.0j)

        I = C + A*numpy.exp((0+1j)*theta) + B*numpy.exp((0+1j)*2.0*theta)
        return numpy.array([I.real, I.imag])

    def errfunc(self, p, theta, spots):
        guessX, guessY = self.fitfunc(p, theta)
        retval = numpy.abs(guessX-spots.T[0]) + numpy.abs(guessY-spots.T[1])
        return retval
        
    def fitLimacons(self, ax=None):
        goodFit = True
        while True:
            #coeffs = random.randn(9).tolist()
            coeffs = random.randn(6).tolist()
            focal_fit = optimize.leastsq(self.errfunc, coeffs, 
                    args=(self.theta, self.focal_spots), factor=0.5)

            #focal_fits = focal_fit[0][0] + focal_fit[0][1]*numpy.cos(self.theta)
            difference = self.errfunc(focal_fit[0], self.theta, self.focal_spots)
            points = self.fitfunc(focal_fit[0], self.theta)
            span = numpy.max(points[0])-numpy.min(points[0])
            if ax!= None:
                ax.scatter(points[0], points[1], color = 'r')
                ax.figure.show()
                print numpy.max(difference)/span
                raw_input()

            if numpy.max(difference)/span < 1e-1:
                break
            goodFit = False
            break

        while True:
            #coeffs = random.randn(9).tolist()
            coeffs = numpy.abs(random.randn(6)).tolist()
            pupil_fit = optimize.leastsq(self.errfunc, coeffs,
                    args=(self.theta, self.pupil_centroids), factor=0.5)
            #pupil_fits = pupil_fit[0][0] + pupil_fit[0][1]*numpy.cos(self.theta)
            difference = self.errfunc(pupil_fit[0], self.theta, self.pupil_centroids)
            points = self.fitfunc(pupil_fit[0], self.theta)
            span = numpy.max(points[0])-numpy.min(points[0])
            if ax!= None:
                ax.scatter(points[0], points[1], color = 'b')
                ax.figure.show()
                dumpFit(pupil_fit[0], 'Pupil')
                raw_input()

            if numpy.max(difference)/span < 1e-1:
                break
            goodFit = False
            break

        if goodFit:
            """
            retval = numpy.append(focal_fit[0], pupil_fit[0])
            retval[1] = numpy.abs(retval[1])
            retval[2] = numpy.abs(retval[2])
            retval[4] = numpy.abs(retval[4])
            retval[5] = numpy.abs(retval[5])
            retval[7] = numpy.abs(retval[7])
            retval[8] = numpy.abs(retval[8])
            retval[10] = numpy.abs(retval[10])
            retval[11] = numpy.abs(retval[11])
            dumpFit(retval[0:6], 'Field')
            dumpFit(retval[6:], 'Pupil')
            raw_input()
            return numpy.array(retval)
            #"""
            #"""
            return numpy.append(numpy.append(numpy.append(numpy.append(numpy.append(
                    numpy.append(numpy.append(focal_fit[0][0], focal_fit[0][2]),
                    focal_fit[0][3]),
                    focal_fit[0][5]), pupil_fit[0][0]), pupil_fit[0][2]),
                    pupil_fit[0][3]),
                    pupil_fit[0][5])
            #"""
        else:
            return numpy.zeros(6)

def dumpFit(fit, name):
    print("%s X: 0 = %.2f, 1 = %.2f, 2 = %.2f" % (name, fit[0], fit[1], fit[2]))
    print("%s Y: 0 = %.2f, 1 = %.2f, 2 = %.2f" % (name, fit[3], fit[4], fit[5]))

class Point( object ):
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def __add__(self, other):
        return Point(self.x+other.x, self.y+other.y, self.z+other.z)

    def __sub__(self, other):
        return Point(self.x-other.x, self.y-other.y, self.z-other.z)

    def __div__(self, factor):
        return Point(self.x/factor, self.y/factor, self.z)

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

        
#pyplot.ion()
fig1 = pyplot.figure(0)
fig1.clear()
fig2 = pyplot.figure(1)
fig2.clear()
ax1 = fig1.add_axes([0.1, 0.1, 0.8, 0.4])
ax2 = fig1.add_axes([0.1, 0.5, 0.8, 0.4])
ax3 = fig2.add_axes([0.1, 0.1, 0.8, 0.8])

dX = 0.5
dAngle = 0.1

modes = []
focus_measurements = []
pupil_measurements = []
measurements = []


for i in range(1000):
    print i
    internalMisalignments = {}
    internalMisalignments["m1Point_Z"] = random.randn()*dX
    internalMisalignments["m2Point_Y"] = random.randn()*dX
    internalMisalignments["m3Point_Z"] = random.randn()*dX
    internalMisalignments["m1Angle"] = random.randn()*dAngle
    internalMisalignments["m1Axis"] = random.rand(3)*2.0 - 1.0
    internalMisalignments["m2Angle"] = random.rand()*dAngle
    internalMisalignments["m2Axis"] = random.rand(3)*2.0 - 1.0
    internalMisalignments["m3Angle"] = random.rand()*dAngle
    internalMisalignments["m3Axis"] = random.rand(3)*2.0 - 1.0
    """
    internalMisalignments["m1Point_Z"] = 0.0
    internalMisalignments["m2Point_Y"] = 0.0
    internalMisalignments["m3Point_Z"] = 0.0
    internalMisalignments["m1Angle"] = 0.0
    internalMisalignments["m1Axis"] = [0.0, 0.0, 0.0]
    internalMisalignments["m2Angle"] = 0.0
    internalMisalignments["m2Axis"] = [0.0, 0.0, 0.0]
    internalMisalignments["m3Angle"] = 0.0
    internalMisalignments["m3Axis"] = [0.0, 0.0, 0.0]
    #"""

    tipAngle = random.randn()*dAngle
    tipAxis = numpy.array([1.0, 0.0, 0.0])
    tiltAngle = random.randn()*dAngle
    tiltAxis = numpy.array([0.0, 1.0, 0.0])
    x = random.rand()*dX
    y = random.rand()*dX
    z = 0.0

    derot = Derotator(iMs=internalMisalignments, tipAngle=tipAngle, tiltAngle=tiltAngle,
            tipAxis = tipAxis, tiltAxis=tiltAxis, x=x, y=y, z=z, ax=None)

    fit = derot.fitLimacons(ax=None)
    if numpy.sum(fit) != 0.0:
        modes.append([tipAngle, tiltAngle, x, y])
        measurements.append(fit)

modes = numpy.array(modes)

measurements = numpy.array(measurements)

IM = measurements.T.dot(scipy.linalg.pinv(modes).T)

U, S, V = scipy.linalg.svd(IM)
CM = scipy.linalg.pinv(IM)
#ax.scatter(derot.spots.T[0], derot.spots.T[1], color = 'r')
#ax.scatter(guessX, guessY, color = 'g')
#guessX, guessY = derot.fitfunc(parameters, derot.theta)
#for spx, spy, gX, gY in zip(derot.spots.T[0], derot.spots.T[1], guessX, guessY):
#    ax.plot([spx, gX], [spy, gY], color = 'k')
#ax.imshow(measurements)
#fig.colorbar(ax.matshow(CM))
#ax.set_aspect('equal')
ax1.matshow(V)
fig1.colorbar(ax2.matshow(U))
ax3.plot(S)
fig1.show()
fig2.show()
        #print retval
#newx = []
#newy = []
#for r, theta in zip(fits, derot.theta):
#    newx.append(r*numpy.cos(theta+3.14159))
#    newy.append(r*numpy.sin(theta+3.14159))

#ax.scatter(newx, newy, color = 'r')
#fig.show()
