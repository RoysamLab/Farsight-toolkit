#!/usr/bin/env python
#############################################################################
# A SCRIPT TO RUN THROUGH THE 5 LABEL TEST IMAGE FROM BEGINNING TO END:
################################################################################
import subprocess
import os

from farsightutils import *

###############################################################################

#SET WOKRING DIRECTORY WHERE FILES WILL BE SAVED
SetWorkingDirectory()

image = 'Montage5Funmixed02.lsm'
input_nuc_image = 'Montage5Funmixed02_Nuclei.tiff'
nuc_label_image = 'Montage5Funmixed02_Nuclei_label.tiff'
parameters = 'Seg_Params.ini'
assocXML = 'Montage5Funmixed02_assoc_def.xml'
path = os.getcwd()

#SPLIT THE LSM CHANNELS
subprocess.call(['lsm_to_tiff.exe', image])

#SEGMENT
subprocess.call(['segment_nuclei.exe', input_nuc_image, nuc_label_image, parameters])
#ASSOCIATIVE MEASUREMENTS
subprocess.call(['compute_associative_measures', assocXML])
#FEATURES
subprocess.call(['compute_nuclei_features', path, input_nuc_image, nuc_label_image])

