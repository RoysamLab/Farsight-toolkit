#!/usr/bin/env python
#############################################################################
# HOW TO GET STARTED:
#
# 1. import this module: import associate
# 2. EXECUTE: >>>associate.run()
################################################################################import subprocess
import subprocess
import threading
import shutil
import sys
import os

###############################################################################################
#START ASSOCIATE FUNCTION
def run():
    
  doNewParams = raw_input('Do you want to create a new associative parameters file(y/n)?:  ')

  if doNewParams == 'n':
      
    #print out filenames of XML in this folder
    files = os.listdir(os.getcwd())
    print '\nXML in Working Directory: '
    for f in files:
      if f.find('.xml') != -1:
        print f
    print '\n'

    filename = raw_input('Please choose XML file with association rules:  ')
    filefull = os.getcwd() + os.sep + filename
    ok = os.path.exists(filefull)
    if not ok:
      print 'file: ' + filename + ' does not exist'
      return

    subprocess.call(['compute_associative_measures', filename])
    
  elif doNewParams == 'y':
  
    numRules = raw_input('Please specify number of associative measures to compute:  ')

    #print out filenames of Images in this folder
    files = os.listdir(os.getcwd())
    print '\nImages in Working Directory: '
    for f in files:
      if f.find('.tif') != -1 or f.find('.tiff') != -1 or f.find('.pic') != -1:
        print f
    print '\n'

    filename = raw_input('Please choose input label image file:  ')
    filefull = os.getcwd() + os.sep + filename
    ok = os.path.exists(filefull)
    if not ok:
      print 'file: ' + filename + ' does not exist'
      return

    outName = os.path.splitext(filename)[0] + '_assoc.xml'
    
    subprocess.call(['compute_associative_measures',filename, outName, numRules])

  print '\n\n\n'


#END SEGMENT.NUCLEI FUNCTION
################################################################################

if __name__ == "__main__":
    run()
