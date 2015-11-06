README

How to measure the Pupil and Focus runouts

0) Save your data in the data directory

    cp data.files /home/fprakt/Data/Derotator/derot_directory/

1) Open a terminal.  Type:

    source ~/.tcshrc   (Only necessary until I figure out how to make this automatic (i.e. I ask Udo))

2) change the hard-coded directories and file name structures in the /home/fprakt/Code/Graffity/src/CentroidFinder/ files

   focus_finder.py
     line 14 : directory = "Path/To/DATAFILES/"
     line 15 : files = glob.glob(directory+'focus*001.TIF')  <- Ensure this follows filename format
   pupil_finder.py
     line 16 : directory= '/home/fprakt/Data/Derotator/derot_04112015/'
     line 17 : files = glob.glob(directory+'pupil*001.TIF') <- Ensure this follows filename format

  (If pupil has moved significantly, you might want to change the X and Y values for the pupil given in lines 20-24 of pupil_finder.py)
  
   Save your changes.

3) start ipython

    ipython

4) run focus_finder.py
    - saves the X and Y positions of focus in Focus_Runout.txt

5) run pupil_finder.py
    - saves the X and Y positions of center of pupil in Pupil_Runout.txt

6) run plot_runout.py
    - uses both output files to both runouts on a single .png file.