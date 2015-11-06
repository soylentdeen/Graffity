Zernike Fun-  Manipulating Zernikes on the DM for Fun and Profit!
=====================================================================

Overview:

   Are you sick and tired of trying to align optical systems without full 
control of the deformable mirror?  Are aberrations like trefoil and astigmatism making your life miserable?  Then try MPIA-CIAO's brand new ZernikeFun!  ZernikeFun allows you to send static patterns specified in Zernike amplitudes to the VLT DM to improve your optical system performance!  It slices! It dices!  All while curing male-pattern baldness** and cleaning leaves from your gutter!

ZernikeFun allows you to specify a sequence of Zernike coefficients (in the Noll ordering convention) in units of RMS microns.  ZernikeFun then uses a patent-pending ZERNIKE-2-DM matrix to compute the actuator offsets relative to the quasistatic pattern, and saves the new actuator positions in a results.fits file.  Copying this new .fits file over to the SPARTA machine and applying it to the mirror has never been easier!!  Order today!

**ZernikeFun has not been approved by the FDA to treat baldness, and this advertisement should not be construed as medical advice.  Side-effects may include extreme flatulence, fits of giggling, and/or delusions of grandeur.


How to Use ZernikeFun:
------------------------------
1) Set up the SPARTA system for use with ZernikeFun
   - log in to the SPARTA gateway machine in the ciaomgr account:
       ssh -X ciaomgr@wci1sgw
       
   - Start the SPARTA panel
       wci1sgw ciaomgr: 100> spcipanei &
       
   - Initialize the system, if it is not already INITialized or ONLINE
       - If MAIN state is UNDEFINED:
           - MENU <Startup/Shutdown> -> Start Config CIAO
       - If MAIN state is NOT-INIT:
           - MENU <State> -> INIT
           
   - Open the CODE 

   
   
2) Open ipython and prepare to run ZernikeFun
   - Graffity/ZernikeFun> ipython3 --pylab=tk
   
   
   https://en.wikipedia.org/wiki/Zernike_polynomials