#############################################################################
# A SCRIPT TO RUN THROUGH THE 5 LABEL TEST IMAGE FROM BEGINNING TO END:
################################################################################import subprocess
import subprocess
import threading
import shutil
import sys
import os

###############################################################################################
inputimage = 'Montage5Funmixed02_Nuclei.tiff'
segimage = 'Montage5Funmixed02_Nuclei_label.tiff'
parameters = 'Seg_Params.ini'
path = os.getcwd()

#SEGMENT
subprocess.call(['segment_nuclei.exe', inputimage, segimage, parameters])
#FEATURES
subprocess.call(['compute_nuclei_features', path, inputimage, segimage])


