import scipy
import numpy
import matplotlib.pyplot as pyplot
import pyfits
import py_compile
import sys
import paramiko
import scp
import time
import glob

def createSSHClient(server, port, user, password):
    client = paramiko.SSHClient()
    client.load_system_host_keys()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(server, port, user, password)
    return client

def parseQSP(pattern):
    nAct = int(len(pattern)/4)
    actPos = []
    for i in range(nAct):    
        actPos.append(2.0*(int(pattern[i*4:(i+1)*4], 16)+1.0)/2**14 - 1.0)
    return numpy.array(actPos)

def getChoice():
    while True:
        print("")
        print("ZernikeFun! Manipulating Zernikes on the CIAO DM for Fun and Profit!")
        print("Select on the of the following options:")
        print("  -(Z)ero the DM")
        print("  -Save current actuator pattern in .(F)its file")
        print("  -Save current Zernike coefficients in (T)ext file")
        print("  -(L)oad actuator pattern from .fits file")
        print("  -(D)isplay current Zernike coefficients")
        print("  -Enter (A)bsolute Zernike coefficients")
        print("  -Enter (R)elative Zernike coefficient offsets")
        print("  -(M)ove the Tip/Tilt Mount")
        print("  -(Q)uit")
        print("====================================================================")
        choice = input("Your choice:")
        if choice.upper() in ['Z', 'F', 'T', 'L', 'D', 'A', 'R', 'Q']:
            return choice.upper()


def sendPattern(current, ssh):
    scaled = numpy.array((current+1.0)/2.0 *2.0**14 -1.0, dtype=numpy.int16)

    command = ''
    for act in scaled:
        command += "%0.2X" % act
    
    if not(numpy.any(numpy.abs(current) > 0.98)):
        ssh.exec_command('sciapiSendCommand lci1hva QSPATTERN SET %s' % command)
        pass
    else:
        print("Error!  New Mirror Position will induce saturations!")


server = 'wci1sgw'
port = 22
user = 'ciaomgr'
password = 'ciao_2015'
ssh = createSSHClient(server, port, user, password)
scp = scp.SCPClient(ssh.get_transport())

Z2DM = pyfits.getdata('cimdatZernike2DM.fits')
DM2Z = numpy.linalg.pinv(Z2DM)

nominal = pyfits.getdata('QSPATTERN.fits')

noll_index = numpy.arange(len(Z2DM)) + 1

stdin, stdout, stderr = ssh.exec_command('sciapiSendCommand lci1hva QSPATTERN GET')
current = parseQSP(stdout.readlines()[0].strip())

current_coeffs = current[:60].dot(DM2Z)*1e6

actnum = numpy.arange(60)

fig = pyplot.figure(0)
fig.clear()
ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.4])
ax2 = fig.add_axes([0.1, 0.5, 0.8, 0.4])

ax1.plot([0, 60], [-1, -1], color = 'r')
ax1.plot([0, 60], [1, 1], color = 'r')
nom_plot, = ax1.plot(actnum, nominal[0], label = 'Nominal Flat')
cur_plot, = ax1.plot(actnum, current[:60], label = 'Current')
ax1.set_ybound(lower=-1.1, upper=1.1)

coeff_plot = ax2.bar(noll_index-0.5, current_coeffs)


fig.show()


choice = getChoice()

Tip = 0.0
Tilt = 0.0

while choice != 'Q':
    if choice == 'R':
        print("Enter space-separated pairs of deltas to the Noll Coeffecients (rms microns):")
        print("Example: \'1 0.5 4 0.1\' will add 0.5 um to Piston and 0.1 um to Focus")
        input_string = input("Enter relative offsets: ")
        if input_string == 'Q':
            choice = getChoice()
        else:
            try:
                coeffs = numpy.zeros(len(Z2DM))
                input_string = numpy.array(input_string.split(), dtype=numpy.float32)
                input_indices = input_string[::2]
                input_coeffs = input_string[1::2]
                coeffs[numpy.array(input_indices, dtype=int)-1] = input_coeffs
                delta = numpy.dot(coeffs*1e-6, Z2DM)
                current = current[:60].copy() + delta
    
                cur_plot.set_data(actnum, current)
                ax2.clear()
                ax2.bar(noll_index-0.5, current.dot(DM2Z)*1e6)
    
                fig.canvas.draw()

                # Set TT to Zero
                current = numpy.append(current, Tip)
                current = numpy.append(current, Tilt)
            
                sendPattern(current, ssh)
            except:
                print("Error!")
    elif choice == 'A':
        print("Enter space-separated pairs of Noll Coeffecients (rms microns):")
        print("Example: \'1 0.5 4 0.1\' will put 0.5 um to Piston and 0.1 um to Focus")
        input_string = input("Enter Absolute coefficients: ")
        if input_string == 'Q':
            choice = getChoice()
        else:
            try:
                coeffs = numpy.zeros(len(Z2DM))
                input_string = numpy.array(input_string.split(), dtype=numpy.float32)
                input_indices = input_string[::2]
                input_coeffs = input_string[1::2]
                coeffs[numpy.array(input_indices, dtype=int)-1] = input_coeffs
                current = numpy.dot(coeffs*1e-6, Z2DM)
                
                cur_plot.set_data(actnum, current)
                ax2.clear()
                ax2.bar(noll_index-0.5, current.dot(DM2Z)*1e6)
    
                fig.canvas.draw()
    
                # Set TT to Zero
                current = numpy.append(current, Tip)
                current = numpy.append(current, Tilt)
                
                sendPattern(current, ssh)
            except:
                print("Error!")
    elif choice == 'Z':
        coeffs = numpy.zeros(len(Z2DM))
        current = numpy.zeros(62)
        
        cur_plot.set_data(actnum, current[:60])
        ax2.clear()
        ax2.bar(noll_index-0.5, current[:60].dot(DM2Z)*1e6)

        fig.canvas.draw()
        
        sendPattern(current, ssh)
        choice = getChoice()
    elif choice == 'F':
        print("Here are the current .fits files in the current directory")
        print(glob.glob('*.fits'))
        print("Enter .fits file name to save actuator positions to:")
        filename = input("(Default = ZernFun_Results.fits) :")
        if filename == '':
            filename = "ZernFun_Results.fits"
        HDU = pyfits.PrimaryHDU(current[:60])
        HDU.header.set("DATE", time.ctime())
        for i, z in zip(noll_index, current[:60].dot(DM2Z)*1e6):
            HDU.header.set("NOLL_"+str(i), z)
        HDU.writeto(filename, clobber=True)
        choice = getChoice()
    elif choice == 'T':
        print("Here are the current .txt files in the current directory")
        print(glob.glob('*.txt'))
        print("Enter .txt file name to save Zernike coefficients to:")
        filename = input("(Default = ZernFun_Coefficients.txt) :")
        if filename == '':
            filename = "ZernFun_Coefficients.txt"
        outfile = open(filename, 'w')
        outfile.write("Zernike Coefficients : %s\n" % time.ctime())
        outfile.write("========================================\n")
        for i, z in zip(noll_index, current[:60].dot(DM2Z)*1e6):
            outfile.write("%-15s %-15s\n" % (str(i), '%.3f'%z))
        outfile.close()
        choice = getChoice()
    elif choice == 'L':
        print("Here are the current .fits files in the current directory")
        print(glob.glob('*.fits'))
        filename = input("Enter the name of the .fits from which to read the actuator positions:")
        current = pyfits.getdata(filename)
        current = numpy.append(current, 0.0)
        current = numpy.append(current, 0.0)
        
        cur_plot.set_data(actnum, current[:60])
        ax2.clear()
        ax2.bar(noll_index-0.5, current[:60].dot(DM2Z)*1e6)

        fig.canvas.draw()
        
        sendPattern(current, ssh)
        choice = getChoice()
    elif choice == 'D':
        print("Zernike Coefficients : %s" % time.ctime())
        print("Noll Index     | RMS (microns)")
        print("========================================")
        for i, z in zip(noll_index, current[:60].dot(DM2Z)*1e6):
            print("%-15s | %-15s" % (str(i), '%.3f'%z))
        choice = getChoice()

