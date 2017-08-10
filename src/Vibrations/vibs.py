import sys

sys.path.append('../')

import Graffity
from matplotlib import pyplot
import CIAO_DatabaseTools

CIAO_keywords = ['STREHL', 'SEEING', 'ALT', 'AZ']

CIAO_DB = CIAO_DatabaseTools.CIAO_Database()

CIAO_values = CIAO_DB.query(keywords=CIAO_keywords, timeOfDay='NIGHT', startTime='2017-08-09 00:00:00')

AVCModes = {}

for CIAO_ID in [1]:
    for record in CIAO_values[CIAO_ID]:
        dl = Graffity.DataLogger(directory=record[-3])
        dl.loadData()
        dl.computeStrehl()
        dl.measureVibs(frequencies=[48.0, 99.0, 115.0])
        print asdf



