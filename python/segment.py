#############################################################################
# HOW TO GET STARTED:
#
# 1. import this module: import segment
# 2. SEGMENT: >>>segment.nuclei()
#    YOU WILL BE ASKED TO PROVIDE PATH CONTAINING THE IMAGE AND PARAMETER FILE
#    MAKES SURE PARAMETERS ARE SET BEFORE CALLING THIS FUNCTION
#    THE FUNCTION WILL CREATE A NUMER OF OUTPUT FILES IN THIS SAME FOLDER
################################################################################import subprocess
import subprocess
import threading
import shutil
import sys
import os

###############################################################################################
#START SEGMENT.NUCLEI FUNCTION
def nuclei():
  #print out filenames of images in this folder
  files = os.listdir(os.getcwd())
  print '\nImages in Working Directory are: '
  for f in files:
    if f.find('.tif') != -1 or f.find('.pic') != -1:
      print f
  print '\n'

  filename = raw_input('Please choose input file:  ')

  filefull = os.getcwd() + os.sep + filename
  ok = os.path.exists(filefull)
  if not ok:
    print 'file: ' + filename + ' does not exist'
    return

  #print out .ini files to find parameter files
  print '\nConfiguration Files in Working Directory are: '
  for f in files:
    if f.find('.ini') != -1:
      print f
  print '\n'

  params = raw_input('Please type filename of parameter file(Seg_Params.ini):  ')
  if params == '':
    params = 'Seg_Params.ini';

  parafull = os.getcwd() + os.sep + params
  ok = os.path.exists(parafull)
  if not ok:
    print 'params: ' + params + ' does not exist'
    return
  
  print '\n\n\n'
  
  #create output filename:
  (begin,end) = os.path.splitext(filefull)
  outfull = begin + '_label' + end;
  subprocess.call(['segment_nuclei.exe', filefull, outfull, parafull])

  #NOW COMPUTE THE FEATURES AND CREATE XML
  (mypath,infile) = os.path.split(filefull)
  (mypath,lblfile) = os.path.split(outfull)
  subprocess.call(['compute_nuclei_features.exe',mypath,infile,lblfile])

#END SEGMENT.NUCLEI FUNCTION
###############################################################################################
