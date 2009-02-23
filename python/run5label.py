#!/usr/bin/env python
#############################################################################
# A SCRIPT TO RUN THROUGH THE 5 LABEL TEST IMAGE FROM BEGINNING TO END:
################################################################################
import subprocess
import os

from farsightutils import *

###############################################################################
parameters = 'Seg_Params.ini'
path = os.getcwd()

assocXML = 'w2less_assoc.xml'
numRules = '1'

#GET THE FILENAME OF THE LSM IMAGE
image = GetFilename()

#SET WOKRING DIRECTORY WHERE FILES WILL BE SAVED
SetWorkiingDirectory()

#SPLIT THE LSM CHANNELS
subprocess.call(['lsm_to_tiff.exe', image])

#SEGMENT
#subprocess.call(['segment_nuclei.exe', inputimage, segimage, parameters])
#FEATURES
#subprocess.call(['compute_nuclei_features', path, inputimage, segimage])
#ASSOCIATIVE MEASUREMENTS
#subprocess.call(['compute_associative_measures',segimage, assocXML, numRules])

