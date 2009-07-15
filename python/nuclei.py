#############################################################################
# HOW TO GET STARTED:
#
# 1. OPEN THIS FILE (nuclei.py) USING IDLE
# 2. IMPORT THIS MODULE: >>>import nuclei
# 3. SEGMENT: >>>nuclei.segment()
#    YOU WILL BE ASKED TO PROVIDE PATH CONTAINING THE IMAGE AND PARAMETER FILE
#    MAKES SURE PARAMETERS ARE SET BEFORE CALLING THIS FUNCTION
#    THE FUNCTION WILL CREATE A NUMER OF OUTPUT FILES IN THIS SAME FOLDER
# 4. OUTLIER DETECTOR:  >>>nuclei.outlier()
#    INPUT FILE IS CALLED imgname_features.txt
# 5. START GUI:  >>>nuclei.show()
# 6. LOAD RESULT: IN GUI File->Open Result
#     CHANGE 'COLOR BY' FIELD IN PLOT TO SEE OUTLIERS
################################################################################import subprocess
import subprocess
import threading
import shutil
import sys
import os

###############################################################################################
#START SEGMENT FUNCTION
def segment():
  inpath = raw_input('Please specify path to working directory:  ')
  if not os.path.exists(inpath):
    print 'Path: ' + inpath + ' does not exist'
    return

  #print out filenames of images in this folder
  files = os.listdir(inpath)
  print '\n'
  for f in files:
    if f.find('.tif') != -1 or f.find('.pic') != -1:
      print f
  print '\n'

  filename = raw_input('Please type filename from list above:  ')

  filefull = inpath + os.sep + filename
  ok = os.path.exists(filefull)
  if not ok:
    print 'file: ' + filename + ' does not exist'
    return

  #print out .ini files to find parameter files
  print '\n'
  for f in files:
    if f.find('.ini') != -1:
      print f
  print '\n'

  params = raw_input('Please type filename of parameter file(Seg_Params.ini):  ')
  if params == '':
    params = 'Seg_Params.ini';

  parafull = inpath + os.sep + params
  ok = os.path.exists(parafull)
  if not ok:
    print 'params: ' + params + ' does not exist'
    return
  
  print '\n\n\n'
  subprocess.call(['nucseg2.exe', filefull, parafull])

#END SEGMENT FUNCTION
###############################################################################################
###############################################################################################
#START SHOW FUNCTION
class FarsightThread ( threading.Thread ):
   def run ( self ):
      subprocess.call(["Farsight.exe"])

def show():
  print 'Starting Farsight GUI'
  FarsightThread().start()
#END SHOW FUNCTION
###############################################################################################
###############################################################################################
#START OUTLIER FUNCTION
def outlier():
  filename = raw_input('Specify filename of file in metaNeural format:  ')
  if not os.path.exists(filename):
    print 'file: ' + filename + ' does not exist'

  fdir = os.path.abspath( os.path.dirname(filename) ) #absolute directory of the file
  fnam = os.path.basename(filename)			 #name of file

  wrkdir = fdir + os.sep + 'temp'
  if os.path.exists(wrkdir):		#delete all of the files in the directory
    files = os.listdir(wrkdir)
    for f in files:
     	os.remove(wrkdir + os.sep + f)
  else:
    os.mkdir(wrkdir)			#makes a work directory to put all files into
  
  newname = 'cell.txt'
  shutil.copy(filename, wrkdir + os.sep + newname )		#copy the datafile into the work directory & rename

  #change my working directory
  origdir = os.getcwd()
  os.chdir(wrkdir);

  print "\nOUTLIER DETECTION BEGINNING..."
  cout = file('cout.txt','w')
  cin = file('cin.txt','w+')
  cin.write("10\n")
  cin.write("0.15053\n")
  #SCALE DATA AND RUN ONE-CLASS SVM:
  subprocess.Popen(['dmak', newname, '-3'],stdout=cout).wait()
  shutil.copy( 'la_sscala.txt', 'la_ssscala.txt' )
  subprocess.Popen(['analyze', newname + '.txt', '6044'], stdout=cout, stdin=cin).wait()
  cout.close()
  cin.close()
  print "...DONE\n"

  #copy the outlier file back to working directory
  shutil.copy( 'outliers.txt', fdir + os.sep + 'outliers.txt' )
  
  #change back to original working directory
  os.chdir(origdir)

  #print "\nOutlier list has been saved to project directory"
  #print "Type nuclei.show() to load Farsight Viewer"
  #print "Then Load the XML file to view the results\n"
#END OUTLIER FUNCTION
################################################################################################
