#############################################################################
# A SCRIPT TO RUN THROUGH THE 5 LABEL TEST IMAGE FROM BEGINNING TO END:
################################################################################
import subprocess
import threading
import shutil
import sys
import os

###############################################################################
inputimage = 'w2less.tif'
segimage = 'w2less_label.tif'
parameters = 'Seg_Params.ini'
path = os.getcwd()

assocXML = 'w2less_assoc.xml'
numRules = '1'

#SEGMENT
subprocess.call(['segment_nuclei.exe', inputimage, segimage, parameters])
#FEATURES
subprocess.call(['compute_nuclei_features', path, inputimage, segimage])
#ASSOCIATIVE MEASUREMENTS
subprocess.call(['compute_associative_measures',segimage, assocXML, numRules])

