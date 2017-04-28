import scipy
import numpy
from matplotlib import pyplot
import Graffity
import glob
import time
import sqlite3

connection = sqlite3.connect('/home/cdeen/Data/CIAO/Dataloggers.db')

cursor = connection.cursor()

fig0 = pyplot.figure(0)
fig0.clear()
ax0 = fig0.add_axes([0.1, 0.1, 0.8, 0.8])

CIAO_dir = '/home/cdeen/Data/CIAO/JanComm/2017-01-08_'#/CIAO1/DATA_LOGGER-000038/'
#CIAO_dir = '/home/cdeen/Data/CIAO/2016-09-21/CIAO'#/CIAO1/DATA_LOGGER-000038/'
CDMS_Dir = '/home/cdeen/Code/CIAO/SPARTA/SPARTA_CIAO/CONFIG/spcicfg/config/RTCDATA/CIAO/DEFAULT/'
CDMS_BaseDir = '/home/cdeen/Code/CIAO/SPARTA/SPARTA_CIAO/CONFIG/spcicfg'
CDMS_ConfigDir = '/config/RTCDATA/CIAO/DEFAULT/'
loopRate = 500 #Hz
RTC_Delay = 0.5e-3

Seeing = {}
ASM_Seeing = {}
Strehl = {}
RMS = {}
Tau0 = {}
Time = {}
Alt = {}
Az = {}
Flux = {}
TTresid = {}
#Parallactic = {}
#Isoplanactic = {}


startTime = time.mktime(time.strptime('2017-01-09T00:57:00', '%Y-%m-%dT%H:%M:%S'))


for ciao in [1,2,3,4]:
    CIAO_datadir = CIAO_dir+str(ciao)+'/'
    seeing = []
    asm_seeing = []
    strehl = []
    rms = []
    tau0 = []
    t = []
    alt = []
    az = []
    flux = []
    ttr = []
    for dp, dn, fn in glob.os.walk(CIAO_datadir):
        if (len(dn) == 0) and ('DATA_LOGGER' in dp):
            print dp
            CB =Graffity.CircularBuffer(dp+'/CIAO_LOOP_0001.fits', CDMS_BaseDir=CDMS_BaseDir, 
                 CDMS_ConfigDir = CDMS_ConfigDir, S2M=dp+'/RecnOptimiser.S2M_0001.fits', 
                 ModalBasis='RecnOptimiser.ModalBasis.fits',Z2DM='RecnOptimiser.Z2DM.fits', 
                 CM=dp+'/Recn.REC1.CM_0001.fits', HOIM=dp+'/RecnOptimiser.HO_IM_0001.fits',
                 TT2HO='RecnOptimiser.TT2HO.fits', DM2Z='RecnOptimiser.DM2Z.fits',
                 TTM2Z='RecnOptimiser.TTM2Z.fits', loopRate=loopRate, RTC_Delay=RTC_Delay)

            
            tm = (time.mktime(time.strptime(CB.header.get('ESO TPL START'), '%Y-%m-%dT%H:%M:%S'))-startTime)/60.0
            currentTime = tm*60.0 + startTime 
            #sqlCommand = "PRAGMA index_list(CIAO_1_DataLoggers);"
            #sqlCommand = "PRAGMA index_info(sqlite_autoindex_CIAO_1_DataLoggers_1);"
            sqlCommand = "SELECT * FROM CIAO_%d_DataLoggers WHERE TIMESTAMP = %.1f" % (ciao, currentTime)

            cursor.execute(sqlCommand)
            result = cursor.fetchall()
            if len(result) == 0:
                CB.synchronizeData()
                CB.zernikeSpace()
                CB.computePSD(source = 'ZSlopes')

                CB.computePSD(source = 'ZCommands')

                CB.AORejectionFunction()
                CB.combinePSDs()
                CB.computeKolmogorovCovar()
                CB.zernikePropagation()
                CB.noiseEvaluation()
                CB.seeingEstimation()
                CB.computeSpectralSlope((3,15))
                CB.estimateStrehlRatio()
                CB.computeTTResiduals()


                print("CIAO #%d"% ciao)
                print("%d" % CB.CIAO_ID)
                print("file: %s " % dp)
                print("Seeing : %.3f" % (CB.Seeing*CB.Arcsec))
                print("Strehl : %.3f" % CB.Strehl)
                print("Tau0   : %.4f" % CB.Tau0)
                seeing.append(CB.Seeing*CB.Arcsec)
                asm_seeing.append(CB.header.get('ESO TEL AMBI FWHM'))
                strehl.append(CB.Strehl)
                rms.append(CB.rms)
                tau0.append(CB.Tau0)
                t.append(CB.header.get('ESO TPL START'))
                alt.append(CB.Alt)
                az.append(CB.Az)
                flux.append(numpy.mean(CB.Intensities))
                ttr.append(numpy.std(CB.TTResiduals, axis=0)*0.5)
                obsTime = tm*60+startTime

                
                values = {}
                values['timestamp'] = obsTime
                values['filename'] = CB.df
                values['seeing'] = seeing[-1]
                values['asm_seeing'] = asm_seeing[-1]
                values['strehl'] = strehl[-1]
                values['tau0'] = tau0[-1]
                values['terr'] = CB.TemporalError
                values['avc_state'] = CB.header.get('ESO AOS AVC LOOP ST') == True
                values['cm_modes'] = CB.header.get('ESO AOS CM MODES CONTROLLED')
                values['wfs_geom'] = CB.header.get('ESO AOS GEOM NAME')
                values['gain'] = CB.header.get('ESO AOS GLOBAL GAIN')
                values['awf_enable'] = CB.header.get('ESO AOS HOCTR AWF ENABLE')
                values['garbage_gain'] = CB.header.get('ESO AOS HOCTR GARBAGE GAIN')
                values['ki'] = CB.header.get('ESO AOS HOCTR KI')
                values['kt'] = CB.header.get('ESO AOS HOCTR KT')
                values['pra_enable'] = CB.header.get('ESO AOS HOCTR PRA ENABLE')
                values['pra_gain'] = CB.header.get('ESO AOS HOCTR PRA GAIN')
                values['sma_enable'] = CB.header.get('ESO AOS HOCTR SMA ENABLE')
                values['sma_high'] = CB.header.get('ESO AOS HOCTR SMA HIGH')
                values['sma_iterations'] = CB.header.get('ESO AOS HOCTR SMA ITERATIONS')
                values['sma_low'] = CB.header.get('ESO AOS HOCTR SMA LOW')
                values['tt_ki'] = CB.header.get('ESO AOS HOCTR TT KI')
                values['tt_kt'] = CB.header.get('ESO AOS HOCTR TT KT')
                values['looprate'] = CB.header.get('ESO AOS LOOP RATE')
                values['tt_loop_state'] = CB.header.get('ESO AOS TT LOOP ST')
                values['ttx_refpos'] = CB.header.get('ESO AOS TTM REFPOS X')
                values['tty_refpos'] = CB.header.get('ESO AOS TTM REFPOS Y')
                values['vib_sr'] = CB.header.get('ESO AOS VIB SR')
                values['wfs_mode'] = CB.header.get('ESO AOS WFS MODE')
                values['im_trk_mode'] = CB.header.get('ESO OCS IM TRK MODE')
                values['alt'] = CB.header.get('ESO TEL ALT')
                values['az'] = CB.header.get('ESO TEL AZ')
                values['dit'] = CB.header.get('ESO DET DIT')
                values['derot_enc'] = CB.header.get('ESO INS DROT ENC')
                values['filt_enc'] = CB.header.get('ESO INS FILT1 ENC')
                values['msel_enc'] = CB.header.get('ESO INS MSEL ENC')
                values['msel_name'] = CB.header.get('ESO INS MSEL NAME')
                values['pmtil_enc'] = CB.header.get('ESO INS PMTIL ENC')
                values['pmtip_enc'] = CB.header.get('ESO INS PMTIP ENC')
                values['fldlx'] = CB.header.get('ESO INS FLDL X')
                values['fldly'] = CB.header.get('ESO INS FLDL Y')
                values['fsm_a_u'] = CB.header.get('ESO STS FSM1 GUIDE U')
                values['fsm_a_w'] = CB.header.get('ESO STS FSM1 GUIDE W')
                values['fsm_a_x'] = CB.header.get('ESO STS FSM1 POSX')
                values['fsm_a_y'] = CB.header.get('ESO STS FSM1 POSY')
                values['fsm_b_u'] = CB.header.get('ESO STS FSM2 GUIDE U')
                values['fsm_b_w'] = CB.header.get('ESO STS FSM2 GUIDE W')
                values['fsm_b_x'] = CB.header.get('ESO STS FSM2 POSX')
                values['fsm_b_y'] = CB.header.get('ESO STS FSM2 POSY')
                values['vcm_a_u'] = CB.header.get('ESO STS VCM1 GUIDE U')
                values['vcm_a_w'] = CB.header.get('ESO STS VCM1 GUIDE W')
                values['vcm_a_x'] = CB.header.get('ESO STS VCM1 POSX')
                values['vcm_a_y'] = CB.header.get('ESO STS VCM1 POSY')
                values['vcm_b_u'] = CB.header.get('ESO STS VCM2 GUIDE U')
                values['vcm_b_w'] = CB.header.get('ESO STS VCM2 GUIDE W')
                values['vcm_b_x'] = CB.header.get('ESO STS VCM2 POSX')
                values['vcm_b_y'] = CB.header.get('ESO STS VCM2 POSY')
                values['m10'] = CB.header.get('ESO STS M10 POSANG')
                values['rhum'] = CB.header.get('ESO TEL AMBI RHUM')
                values['temp'] = CB.header.get('ESO TEL AMBI TEMP')
                values['theta0'] = CB.header.get('ESO TEL AMBI THETA0')
                values['winddir'] = CB.header.get('ESO TEL AMBI WINDDIR')
                values['windsp'] = CB.header.get('ESO TEL AMBI WINDSP')
                values['prltic'] = CB.header.get('ESO TEL PRLTIC')
                values['pup_trk'] = CB.header.get('ESO OCS PUP TRK ST')
                values['flux'] = flux[-1]
                values['tip'] = ttr[-1][0]
                values['tilt'] = ttr[-1][1]

                
                if ciao == 1:
                    format_str = """INSERT INTO CIAO_1_DataLoggers (TIMESTAMP, FILENAME, SEEING, 
                             ASM_SEEING, STREHL, TAU0, TERR, AVC_STATE, CM_MODES, WFS_GEOM, GAIN,
                             HOCTR_AWF_ENABLE, HOCTR_GARBAGE_GAIN, HOCTR_KI, HOCTR_KT, 
                             HOCTR_PRA_ENABLE, HOCTR_PRA_GAIN, HOCTR_SMA_ENABLE, HOCTR_SMA_HIGH,
                             HOCTR_SMA_ITERATIONS, HOCTR_SMA_LOW, HOCTR_TT_KI, HOCTR_TT_KT,
                             LOOPRATE, TT_LOOP_STATE, TTX_REFPOS, TTY_REFPOS, VIB_SR, WFS_MODE,
                             IM_TRK_MODE, ALT, AZ, DIT, DROT_ENC, FILT_ENC, MSEL_ENC, MSEL_NAME,
                             PMTIL_ENC, PMTIP_ENC, FLDLX, FLDLY, FSM_A_U, FSM_A_W, FSM_A_X, FSM_A_Y,
                             FSM_B_U, FSM_B_W, FSM_B_X, FSM_B_Y, VCM_A_U, VCM_A_W, VCM_A_X, VCM_A_Y,
                             VCM_B_U, VCM_B_W, VCM_B_X, VCM_B_Y, M10_POSANG, RHUM, TEMP, THETA0, 
                             WINDDIR, WINDSP, PRLTIC, PUP_TRK, FLUX, TIP_RESIDUALS, TILT_RESIDUALS)
                             VALUES ('{timestamp}', '{filename}', '{seeing}', '{asm_seeing}', 
                             '{strehl}', '{tau0}', '{terr}', '{avc_state}', '{cm_modes}', 
                             '{wfs_geom}', '{gain}', '{awf_enable}', '{garbage_gain}', '{ki}',
                             '{kt}', '{pra_enable}', '{pra_gain}', '{sma_enable}', '{sma_high}',
                             '{sma_iterations}', '{sma_low}', '{tt_ki}', '{tt_kt}', '{looprate}',
                             '{tt_loop_state}', '{ttx_refpos}', '{tty_refpos}', '{vib_sr}',
                             '{wfs_mode}', '{im_trk_mode}', '{alt}', '{az}', '{dit}', '{derot_enc}',
                             '{filt_enc}', '{msel_enc}', '{msel_name}', '{pmtil_enc}', '{pmtip_enc}',
                             '{fldlx}', '{fldly}', '{fsm_a_u}', '{fsm_a_w}', '{fsm_a_x}', '{fsm_a_y}',
                             '{fsm_b_u}', '{fsm_b_w}', '{fsm_b_x}', '{fsm_b_y}', '{vcm_a_u}',
                             '{vcm_a_w}', '{vcm_a_x}', '{vcm_a_y}', '{vcm_b_u}', '{vcm_b_w}',
                             '{vcm_b_x}', '{vcm_b_y}', '{m10}', '{rhum}', '{temp}', '{theta0}',
                             '{winddir}', '{windsp}', '{prltic}', '{pup_trk}', '{flux}', '{tip}',
                             '{tilt}');"""
                elif ciao == 2:
                    format_str = """INSERT INTO CIAO_2_DataLoggers (TIMESTAMP, FILENAME, SEEING, 
                             ASM_SEEING, STREHL, TAU0, TERR, AVC_STATE, CM_MODES, WFS_GEOM, GAIN,
                             HOCTR_AWF_ENABLE, HOCTR_GARBAGE_GAIN, HOCTR_KI, HOCTR_KT, 
                             HOCTR_PRA_ENABLE, HOCTR_PRA_GAIN, HOCTR_SMA_ENABLE, HOCTR_SMA_HIGH,
                             HOCTR_SMA_ITERATIONS, HOCTR_SMA_LOW, HOCTR_TT_KI, HOCTR_TT_KT,
                             LOOPRATE, TT_LOOP_STATE, TTX_REFPOS, TTY_REFPOS, VIB_SR, WFS_MODE,
                             IM_TRK_MODE, ALT, AZ, DIT, DROT_ENC, FILT_ENC, MSEL_ENC, MSEL_NAME,
                             PMTIL_ENC, PMTIP_ENC, FLDLX, FLDLY, FSM_A_U, FSM_A_W, FSM_A_X, FSM_A_Y,
                             FSM_B_U, FSM_B_W, FSM_B_X, FSM_B_Y, VCM_A_U, VCM_A_W, VCM_A_X, VCM_A_Y,
                             VCM_B_U, VCM_B_W, VCM_B_X, VCM_B_Y, M10_POSANG, RHUM, TEMP, THETA0, 
                             WINDDIR, WINDSP, PRLTIC, PUP_TRK, FLUX, TIP_RESIDUALS, TILT_RESIDUALS)
                             VALUES ('{timestamp}', '{filename}', '{seeing}', '{asm_seeing}', 
                             '{strehl}', '{tau0}', '{terr}', '{avc_state}', '{cm_modes}', 
                             '{wfs_geom}', '{gain}', '{awf_enable}', '{garbage_gain}', '{ki}',
                             '{kt}', '{pra_enable}', '{pra_gain}', '{sma_enable}', '{sma_high}',
                             '{sma_iterations}', '{sma_low}', '{tt_ki}', '{tt_kt}', '{looprate}',
                             '{tt_loop_state}', '{ttx_refpos}', '{tty_refpos}', '{vib_sr}',
                             '{wfs_mode}', '{im_trk_mode}', '{alt}', '{az}', '{dit}', '{derot_enc}',
                             '{filt_enc}', '{msel_enc}', '{msel_name}', '{pmtil_enc}', '{pmtip_enc}',
                             '{fldlx}', '{fldly}', '{fsm_a_u}', '{fsm_a_w}', '{fsm_a_x}', '{fsm_a_y}',
                             '{fsm_b_u}', '{fsm_b_w}', '{fsm_b_x}', '{fsm_b_y}', '{vcm_a_u}',
                             '{vcm_a_w}', '{vcm_a_x}', '{vcm_a_y}', '{vcm_b_u}', '{vcm_b_w}',
                             '{vcm_b_x}', '{vcm_b_y}', '{m10}', '{rhum}', '{temp}', '{theta0}',
                             '{winddir}', '{windsp}', '{prltic}', '{pup_trk}', '{flux}', '{tip}',
                             '{tilt}');"""
                elif ciao == 3:
                    format_str = """INSERT INTO CIAO_3_DataLoggers (TIMESTAMP, FILENAME, SEEING, 
                             ASM_SEEING, STREHL, TAU0, TERR, AVC_STATE, CM_MODES, WFS_GEOM, GAIN,
                             HOCTR_AWF_ENABLE, HOCTR_GARBAGE_GAIN, HOCTR_KI, HOCTR_KT, 
                             HOCTR_PRA_ENABLE, HOCTR_PRA_GAIN, HOCTR_SMA_ENABLE, HOCTR_SMA_HIGH,
                             HOCTR_SMA_ITERATIONS, HOCTR_SMA_LOW, HOCTR_TT_KI, HOCTR_TT_KT,
                             LOOPRATE, TT_LOOP_STATE, TTX_REFPOS, TTY_REFPOS, VIB_SR, WFS_MODE,
                             IM_TRK_MODE, ALT, AZ, DIT, DROT_ENC, FILT_ENC, MSEL_ENC, MSEL_NAME,
                             PMTIL_ENC, PMTIP_ENC, FLDLX, FLDLY, FSM_A_U, FSM_A_W, FSM_A_X, FSM_A_Y,
                             FSM_B_U, FSM_B_W, FSM_B_X, FSM_B_Y, VCM_A_U, VCM_A_W, VCM_A_X, VCM_A_Y,
                             VCM_B_U, VCM_B_W, VCM_B_X, VCM_B_Y, M10_POSANG, RHUM, TEMP, THETA0, 
                             WINDDIR, WINDSP, PRLTIC, PUP_TRK, FLUX, TIP_RESIDUALS, TILT_RESIDUALS)
                             VALUES ('{timestamp}', '{filename}', '{seeing}', '{asm_seeing}', 
                             '{strehl}', '{tau0}', '{terr}', '{avc_state}', '{cm_modes}', 
                             '{wfs_geom}', '{gain}', '{awf_enable}', '{garbage_gain}', '{ki}',
                             '{kt}', '{pra_enable}', '{pra_gain}', '{sma_enable}', '{sma_high}',
                             '{sma_iterations}', '{sma_low}', '{tt_ki}', '{tt_kt}', '{looprate}',
                             '{tt_loop_state}', '{ttx_refpos}', '{tty_refpos}', '{vib_sr}',
                             '{wfs_mode}', '{im_trk_mode}', '{alt}', '{az}', '{dit}', '{derot_enc}',
                             '{filt_enc}', '{msel_enc}', '{msel_name}', '{pmtil_enc}', '{pmtip_enc}',
                             '{fldlx}', '{fldly}', '{fsm_a_u}', '{fsm_a_w}', '{fsm_a_x}', '{fsm_a_y}',
                             '{fsm_b_u}', '{fsm_b_w}', '{fsm_b_x}', '{fsm_b_y}', '{vcm_a_u}',
                             '{vcm_a_w}', '{vcm_a_x}', '{vcm_a_y}', '{vcm_b_u}', '{vcm_b_w}',
                             '{vcm_b_x}', '{vcm_b_y}', '{m10}', '{rhum}', '{temp}', '{theta0}',
                             '{winddir}', '{windsp}', '{prltic}', '{pup_trk}', '{flux}', '{tip}',
                             '{tilt}');"""
                elif ciao == 4:
                    format_str = """INSERT INTO CIAO_4_DataLoggers (TIMESTAMP, FILENAME, SEEING, 
                             ASM_SEEING, STREHL, TAU0, TERR, AVC_STATE, CM_MODES, WFS_GEOM, GAIN,
                             HOCTR_AWF_ENABLE, HOCTR_GARBAGE_GAIN, HOCTR_KI, HOCTR_KT, 
                             HOCTR_PRA_ENABLE, HOCTR_PRA_GAIN, HOCTR_SMA_ENABLE, HOCTR_SMA_HIGH,
                             HOCTR_SMA_ITERATIONS, HOCTR_SMA_LOW, HOCTR_TT_KI, HOCTR_TT_KT,
                             LOOPRATE, TT_LOOP_STATE, TTX_REFPOS, TTY_REFPOS, VIB_SR, WFS_MODE,
                             IM_TRK_MODE, ALT, AZ, DIT, DROT_ENC, FILT_ENC, MSEL_ENC, MSEL_NAME,
                             PMTIL_ENC, PMTIP_ENC, FLDLX, FLDLY, FSM_A_U, FSM_A_W, FSM_A_X, FSM_A_Y,
                             FSM_B_U, FSM_B_W, FSM_B_X, FSM_B_Y, VCM_A_U, VCM_A_W, VCM_A_X, VCM_A_Y,
                             VCM_B_U, VCM_B_W, VCM_B_X, VCM_B_Y, M10_POSANG, RHUM, TEMP, THETA0, 
                             WINDDIR, WINDSP, PRLTIC, PUP_TRK, FLUX, TIP_RESIDUALS, TILT_RESIDUALS)
                             VALUES ('{timestamp}', '{filename}', '{seeing}', '{asm_seeing}', 
                             '{strehl}', '{tau0}', '{terr}', '{avc_state}', '{cm_modes}', 
                             '{wfs_geom}', '{gain}', '{awf_enable}', '{garbage_gain}', '{ki}',
                             '{kt}', '{pra_enable}', '{pra_gain}', '{sma_enable}', '{sma_high}',
                             '{sma_iterations}', '{sma_low}', '{tt_ki}', '{tt_kt}', '{looprate}',
                             '{tt_loop_state}', '{ttx_refpos}', '{tty_refpos}', '{vib_sr}',
                             '{wfs_mode}', '{im_trk_mode}', '{alt}', '{az}', '{dit}', '{derot_enc}',
                             '{filt_enc}', '{msel_enc}', '{msel_name}', '{pmtil_enc}', '{pmtip_enc}',
                             '{fldlx}', '{fldly}', '{fsm_a_u}', '{fsm_a_w}', '{fsm_a_x}', '{fsm_a_y}',
                             '{fsm_b_u}', '{fsm_b_w}', '{fsm_b_x}', '{fsm_b_y}', '{vcm_a_u}',
                             '{vcm_a_w}', '{vcm_a_x}', '{vcm_a_y}', '{vcm_b_u}', '{vcm_b_w}',
                             '{vcm_b_x}', '{vcm_b_y}', '{m10}', '{rhum}', '{temp}', '{theta0}',
                             '{winddir}', '{windsp}', '{prltic}', '{pup_trk}', '{flux}', '{tip}',
                             '{tilt}');"""
                sqlCommand = format_str.format(timestamp=values['timestamp'],filename=values['filename'],
                             seeing=values['seeing'], asm_seeing=values['asm_seeing'], 
                             strehl=values['strehl'], tau0=values['tau0'], terr=values['terr'],
                             avc_state=values['avc_state'], cm_modes=values['cm_modes'],
                             wfs_geom=values['wfs_geom'], gain=values['gain'], 
                             awf_enable=values['awf_enable'], garbage_gain=values['garbage_gain'],
                             ki=values['ki'], kt=values['kt'], pra_enable=values['pra_enable'],
                             pra_gain=values['pra_gain'], sma_enable=values['sma_enable'],
                             sma_high=values['sma_high'], sma_iterations=values['sma_iterations'],
                             sma_low=values['sma_low'], tt_ki=values['tt_ki'], 
                             tt_kt=values['tt_kt'], looprate=values['looprate'],
                             tt_loop_state=values['tt_loop_state'], ttx_refpos=values['ttx_refpos'],
                             tty_refpos=values['tty_refpos'], vib_sr=values['vib_sr'],
                             wfs_mode=values['wfs_mode'], im_trk_mode=values['im_trk_mode'],
                             alt=values['alt'], az=values['az'], dit=values['dit'],
                             derot_enc=values['derot_enc'], filt_enc=values['filt_enc'],
                             msel_enc=values['msel_enc'], msel_name=values['msel_name'],
                             pmtil_enc=values['pmtil_enc'], pmtip_enc=values['pmtip_enc'],
                             fldlx=values['fldlx'], fldly=values['fldly'], fsm_a_u=values['fsm_a_u'],
                             fsm_a_w=values['fsm_a_w'], fsm_a_x=values['fsm_a_x'],
                             fsm_a_y=values['fsm_a_y'], fsm_b_u=values['fsm_b_u'],
                             fsm_b_w=values['fsm_b_w'], fsm_b_x=values['fsm_b_x'],
                             fsm_b_y=values['fsm_b_y'], vcm_a_u=values['vcm_a_u'],
                             vcm_a_w=values['vcm_a_w'], vcm_a_x=values['vcm_a_x'],
                             vcm_a_y=values['vcm_a_y'], vcm_b_u=values['vcm_b_u'],
                             vcm_b_w=values['vcm_b_w'], vcm_b_x=values['vcm_b_x'],
                             vcm_b_y=values['vcm_b_y'], m10=values['m10'], rhum=values['rhum'],
                             temp=values['temp'], theta0=values['theta0'], winddir=values['winddir'],
                             windsp=values['windsp'], prltic=values['prltic'], pup_trk=values['pup_trk'],
                             flux=values['flux'], tip=values['tip'], tilt=values['tilt'])
                             
                cursor.execute(sqlCommand)


    order = numpy.argsort(t)
    Seeing[ciao] = numpy.array(seeing)[order]
    ASM_Seeing[ciao] = numpy.array(asm_seeing)[order]
    Strehl[ciao] = numpy.array(strehl)[order]
    RMS[ciao] = numpy.array(rms)[order]
    Tau0[ciao] = numpy.array(tau0)[order]
    Time[ciao] = numpy.array(t)[order]
    Alt[ciao] = numpy.array(alt)[order]
    Az[ciao] = numpy.array(az)[order]
    Flux[ciao] = numpy.array(flux)[order]
    TTresid[ciao] = numpy.array(ttr)[order]

connection.commit()

connection.close()

for ciao in [1, 2, 3, 4]:
    print "CIAO %d" % ciao
    for t, s, st, alt, az, fl, ttr in zip(Time[ciao], Seeing[ciao], Strehl[ciao],
                                     Alt[ciao], Az[ciao], Flux[ciao], TTresid[ciao]):
        print("%s | %.3f | %.3f | %.2f | %.2f | %.1f | %.2f | %.2f" % 
                    (t, s, st, alt, az, fl, ttr[0]*1000.0, ttr[1]*1000.0))

    ax0.scatter(Seeing[ciao], Strehl[ciao])

fig0.show()

"""
#ax0.clear()

A1 = []
A2 = []
A3 = []
A4 = []
for SR1, SR2, SR3, SR4, seeing1, seeing2, seeing3, seeing4, tau1, tau2, tau3, tau4 in zip(sr_ut1, sr_ut2, sr_ut3, sr_ut4, seeing_ut1, seeing_ut2, seeing_ut3, seeing_ut4, tau_ut1, tau_ut2, tau_ut3, tau_ut4):
    A1.append([1.0, SR1, seeing1, tau1])
    A2.append([1.0, SR2, seeing2, tau2])
    A3.append([1.0, SR3, seeing3, tau3])
    A4.append([1.0, SR4, seeing4, tau4])

A1 = numpy.array(A1)
B1 = numpy.linalg.pinv(A1)
fit1 = B1.dot(psf_ut1)

A2 = numpy.array(A2)
B2 = numpy.linalg.pinv(A2)
fit2 = B2.dot(psf_ut2)

A3 = numpy.array(A3)
B3 = numpy.linalg.pinv(A3)
fit3 = B3.dot(psf_ut3)

A4 = numpy.array(A4)
B4 = numpy.linalg.pinv(A4)
fit4 = B4.dot(psf_ut4)

print fit1
print fit2
print fit3
print fit4
ax0.clear()
ax0.scatter(psf_ut1, sr_ut1*fit1[1] + fit1[0] + seeing_ut1*fit1[2] + tau_ut1*fit1[3], color = 'b')
ax0.scatter(psf_ut2, sr_ut2*fit2[1] + fit2[0] + seeing_ut2*fit2[2] + tau_ut2*fit2[3], color = 'r')
ax0.scatter(psf_ut3, sr_ut3*fit3[1] + fit3[0] + seeing_ut3*fit3[2] + tau_ut3*fit3[3], color = 'g')
ax0.scatter(psf_ut4, sr_ut4*fit4[1] + fit4[0] + seeing_ut4*fit4[2] + tau_ut4*fit4[3], color = 'm')
ax0.plot([0.35, 0.55], [0.35, 0.55])
fig0.show()

raw_input()
ax0.clear()

ax0.errorbar(psf_Time, psf_ut1, psf_dut1, color = 'b', label = 'UT1')
ax0.errorbar(psf_Time, psf_ut2, psf_dut2, color = 'g', label = 'UT2')
ax0.errorbar(psf_Time, psf_ut3, psf_dut3, color = 'r', label = 'UT3')
ax0.errorbar(psf_Time, psf_ut4, psf_dut4, color = 'y', label = 'UT4')

ax0.scatter(Time[1], Strehl[1], color = 'b')
ax0.errorbar(psf_Time, sr_ut1, dsr_ut1, color = 'b')
ax0.scatter(psf_Time, sr_ut1*fit1[1]+fit1[0]+seeing_ut1*fit1[2]+tau_ut1*fit1[3], color = 'b')
ax0.scatter(Time[2], Strehl[2], color = 'g')
ax0.errorbar(psf_Time, sr_ut2, dsr_ut2, color = 'g')
ax0.scatter(psf_Time, sr_ut2*fit2[1]+fit2[0]+seeing_ut2*fit2[2]+tau_ut2*fit2[3], color = 'g')
ax0.scatter(Time[3], Strehl[3], color = 'r')
ax0.errorbar(psf_Time, sr_ut3, dsr_ut3, color = 'r')
ax0.scatter(psf_Time, sr_ut3*fit3[1]+fit3[0]+seeing_ut3*fit3[2]+tau_ut3*fit3[3], color = 'r')
ax0.scatter(Time[4], Strehl[4], color = 'y')
ax0.errorbar(psf_Time, sr_ut4, dsr_ut4, color = 'y')
ax0.scatter(psf_Time, sr_ut4*fit4[1]+fit4[0]+seeing_ut4*fit4[2]+tau_ut4*fit4[3], color = 'y')

ax0.set_xbound(0, 4.0)
ax0.set_ybound(0.0, 0.9)
ax0.set_xlabel("Hours since start of Observing")
ax0.set_ylabel("Computed Strehl Ratio")
pyplot.rcParams["font.size"] = 20.0
ax0.legend(loc = 4)

fig0.show()
#"""
"""
ax0.clear()
#ax0.scatter(Seeing[1], RMS[1], color = 'b')
#ax0.scatter(Seeing[2], RMS[2], color = 'g')
#ax0.scatter(Seeing[3], RMS[3], color = 'r')
#ax0.scatter(Seeing[4], RMS[4], color = 'y')
ax0.scatter(Seeing[1][Seeing[1] > 0.2], Strehl[1][Seeing[1] > 0.2], color = 'b', label='UT1')
ax0.scatter(Seeing[2][Seeing[2] > 0.2], Strehl[2][Seeing[2] > 0.2], color = 'g', label='UT2')
ax0.scatter(Seeing[3][Seeing[3] > 0.2], Strehl[3][Seeing[3] > 0.2], color = 'r', label='UT3')
ax0.scatter(Seeing[4][Seeing[4] > 0.2], Strehl[4][Seeing[4] > 0.2], color = 'y', label='UT4')
ax0.set_xlabel("Estimated Seeing (\")")
ax0.set_ylabel("Estimated Strehl Ratio")
ax0.legend(loc=1, scatterpoints=1)
pyplot.rcParams["font.size"] = 15.0
fig0.show()
fig0.savefig("StrehlvSeeing.png")
#"""
