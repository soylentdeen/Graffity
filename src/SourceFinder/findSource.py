import scipy
import numpy
from numpy import sin, cos, arcsin, deg2rad, rad2deg

class Coordinate( object ):
    def __init__(self, ra, dec, pm=''):
        ra = ra.replace(' ', '')
        dec = dec.replace(' ', '')
        self.pm = pm
        if dec[0] == '+':
            self.sign = 1.0
            self.dsign = ' '
        else:
            self.sign = -1.0
            self.dsign = '-'
        self.rah = ra[:2]
        self.ram = ra[2:4]
        self.ras = ra[4:]
        self.deg = dec[1:3]
        self.dm = dec[3:5]
        self.ds = dec[5:-1]
        self.ra = numpy.int(ra[:2])*15 + numpy.float(ra[2:4])*15.0/60.0 + numpy.float(ra[4:])*15.0/3600.0
        self.dec = self.sign*(numpy.int(dec[1:3])+numpy.float(dec[3:5])/60.0 + numpy.float(dec[5:])/3600.0)

        self.lat = -1.0*(24. + 37.0/60.0 + 38.64/3600.0)
        self.LST = 105.0
        self.ha = self.LST - self.ra
        if self.ha < 0:
            self.ha += 360.0
        self.altitude = rad2deg(arcsin(sin(deg2rad(self.dec)) * sin(deg2rad(self.lat)) + 
                               cos(deg2rad(self.dec))*cos(deg2rad(self.lat))*cos(deg2rad(self.ha))))
        
        print self.altitude


    def __repr__(self):
        return "%s %s %s %s%s %s %s %.3f %s" % (self.rah, self.ram, self.ras, self.dsign, self.deg,
                                        self.dm, self.ds, self.altitude, self.pm)

class Star ( object ):
    def __init__(self, coordinate='', vmag1=30.0, vmag2=30.0, spt1='', spt2='', separation=0.0, pm=''):
        self.coordinate = Coordinate(coordinate[:9], coordinate[9:], pm=pm)
        self.vmag1 = vmag1
        self.vmag2 = vmag2
        self.separation = separation
        self.spt1 = spt1.strip().replace('I','').replace('V','').replace(':','').replace('/','')
        if len(spt2) > 0:
            self.spt2 = spt2.strip().replace('I','').replace('V','').replace(':','').replace('/','')
        else:
            self.spt2 = None

    def calcKmag(self, colors):
        if self.spt1 in colors.keys():
            self.kmag1 = self.vmag1 - colors[self.spt1]
        else:
            self.kmag1 = 30.0
        if self.spt2 != None:
            if self.spt2 in colors.keys():
                self.kmag2 = self.vmag2 - colors[self.spt2]
            else:
                self.kmag2 = 30.0
        else:
            self.kmag2 = 30.0

    def inBox(self, c1, c2, c3, c4):
        if (self.coordinate.ra > c1.ra and self.coordinate.ra < c3.ra):
            if (self.coordinate.dec > c1.dec and self.coordinate.dec < c2.dec):
                return True
            else:
                return False
        else:
            return False


colorsFile = open('starcolors.dat', 'r')

colors = {}

for line in colorsFile.readlines():
    if line[0] != '#':
        l = line.split()
        if l[11][0] != '.':
            colors[l[0][:-1]] = numpy.float(l[11])

df = 'wdsweb_summ2.txt'

data = open(df, 'r')

stars = []
for line in data.readlines():
    if line[0].isdigit():
        try:
            spt = line[70:79].split('+')
            spt1 = spt[0]
            if len(spt) > 1:
                spt2 = spt[1]
            else:
                spt2 = ''
            vmag1 = numpy.float(line[58:63])
            vmag2 = numpy.float(line[64:69])
            coord = line[112:]
            pm = line[80:97]
            separation = (numpy.float(line[46:51]) + numpy.float(line[52:57]))/2.0
            star = Star(coordinate=coord, vmag1=vmag1, vmag2=vmag2, spt1=spt1,
                        spt2=spt2, separation=separation, pm=pm)
            star.calcKmag(colors)
            stars.append(star)
        except:
            pass

ra_min = '020000.0'
ra_max = '120000.0'
dec_min = '-600000.0'
dec_max = '+100000.0'
cutoff = 7.0
sep_min = 7.0
sep_max = 15.0

corner1 = Coordinate(ra_min, dec_min)
corner2 = Coordinate(ra_min, dec_max)
corner3 = Coordinate(ra_max, dec_min)
corner4 = Coordinate(ra_max, dec_max)

possible = []

for star in stars:
    if (star.inBox(corner1, corner2, corner3, corner4) and star.kmag1 < cutoff and
        star.separation > sep_min and star.separation < sep_max and star.kmag2 < 20.0):
        possible.append(star)

#for star in possible:
#    print("Kmag1: %.2f Kmag2: % 6.2f Sep:% 6.2f Coord: %s" % (star.kmag1, star.kmag2, star.separation, star.coordinate))

kmag1 = sorted(possible, key=lambda star: star.kmag1)
kmag2 = sorted(possible, key=lambda star: star.kmag2)
Sep = sorted(possible, key=lambda star: star.separation)

for star in kmag2:
    #print("Kmag1: %.2f Kmag2: % 6.2f Sep:% 6.2f Coord: %s" % (star.kmag1, star.kmag2, star.separation, star.coordinate))
    print("| %s | % 6.2f | % 6.2f | % 6.2f |            |" % (star.coordinate, star.kmag1, star.kmag2, star.separation))
    print("+-------------------------+--------+--------+--------+------------+")


