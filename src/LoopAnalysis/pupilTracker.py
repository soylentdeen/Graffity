import openpyxl
import numpy
from matplotlib import pyplot
import time
import datetime

class Pupil( object ):
    def __init__(self):
        self.t = []
        self.x = []
        self.y = []
        self.z = []

    def add(self, d, t, x, y, z):
        self.t.append(time.mktime(datetime.datetime.combine(datetime.datetime(2016,9,d), t).timetuple()))
        self.x.append(x)
        self.y.append(y)
        self.z.append(z)

    def finish(self):
        self.t = numpy.array(self.t)
        self.x = numpy.array(self.x)
        self.y = numpy.array(self.y)
        self.z = numpy.array(self.z)

fig = pyplot.figure(1)
fig.clear()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

df = '/home/cdeen/Data/GRAVITY/2016-09-21/Pupil.2016-09-21.xlsx'

wb = openpyxl.load_workbook(df)

sheet = wb[wb.sheetnames[0]]

UT1 = Pupil()
UT2 = Pupil()
UT3 = Pupil()
UT4 = Pupil()

for row in sheet.rows[1:]:
    UT1.add(row[2].value, row[3].value, row[4].value, row[5].value, row[6].value)
    UT2.add(row[2].value, row[3].value, row[7].value, row[8].value, row[9].value)
    UT3.add(row[2].value, row[3].value, row[10].value, row[11].value, row[12].value)
    UT4.add(row[2].value, row[3].value, row[13].value, row[14].value, row[15].value)

UT1.finish()
UT2.finish()
UT3.finish()
UT4.finish()


ax.scatter((UT1.t[UT1.x > -100] - UT1.t[0])/3600.0, UT1.x[UT1.x > -100], color = 'b')
ax.plot((UT1.t[UT1.y > -100] - UT1.t[0])/3600.0, UT1.y[UT1.y > -100], color = 'g')
ax.plot((UT1.t[UT1.z > -100] - UT1.t[0])/3600.0, UT1.z[UT1.z > -100], color = 'r')
ax.set_xbound(0.0, 4.0)
fig.show()

raw_input()
ax.clear()

ax.scatter((UT2.t[UT2.x > -100] - UT2.t[0])/3600.0, UT2.x[UT2.x > -100], color = 'b')
ax.plot((UT1.t[UT2.y > -100] - UT2.t[0])/3600.0, UT2.y[UT2.y > -100], color = 'g')
ax.plot((UT1.t[UT2.z > -100] - UT2.t[0])/3600.0, UT2.z[UT2.z > -100], color = 'r')
ax.set_xbound(0.0, 4.0)
fig.show()

raw_input()
ax.clear()

ax.scatter((UT3.t[UT3.x > -100] - UT3.t[0])/3600.0, UT3.x[UT3.x > -100], color = 'b')
ax.plot((UT3.t[UT3.y > -100] - UT3.t[0])/3600.0, UT3.y[UT3.y > -100], color = 'g')
ax.plot((UT3.t[UT3.z > -100] - UT3.t[0])/3600.0, UT3.z[UT3.z > -100], color = 'r')
ax.set_xbound(0.0, 4.0)
fig.show()

raw_input()
ax.clear()

ax.scatter((UT4.t[UT4.x > -100] - UT4.t[0])/3600.0, UT4.x[UT4.x > -100], color = 'b')
ax.plot((UT4.t[UT4.y > -100] - UT4.t[0])/3600.0, UT4.y[UT4.y > -100], color = 'g')
ax.plot((UT4.t[UT4.z > -100] - UT4.t[0])/3600.0, UT4.z[UT4.z > -100], color = 'r')
ax.set_xbound(0.0, 4.0)
fig.show()
