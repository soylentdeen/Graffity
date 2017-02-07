import scipy.io as sio

datafile = '/home/cdeen/Data/CIAO/UT4/BadDatabase/UT4/DMModel/DMModel.mat'

MatlabFile = sio.loadmat(datafile, struct_as_record=False, squeeze_me=True)

DM = MatlabFile['DM']


