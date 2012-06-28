#!/usr/bin/env python
# -*- coding: utf-8 -*-

import shutil
import fnmatch
import os
import subprocess
import os.path

MAIN_PATH = '/FSdata/data/DARPA_MOSAICS'
DATA_FOLDER = ['/0113','/0117','/0120','/0123','/0128','/0131','/0323','/0405','/0409','/0410','/0412','/1206']

MACRO_PATH = ''

print "# Create Folder"

for fileName in DATA_FOLDER:
	PATH = MAIN_PATH+fileName
	if not os.path.isdir(PATH+'_NRRD'):
		print 'creating folder: '+PATH+'_NRRD'
		os.makedirs(PATH+'_NRRD')
	else:
		print 'erasing folder: '+PATH+'_NRRD'
		#for the_file in os.listdir(LOCAL_DATASET_PATH_SEGM_DEBUG):
			#file_path = os.path.join(LOCAL_DATASET_PATH_SEGM_DEBUG, the_file)
			#try:
				#os.unlink(file_path)
			#except Exception, e:
				#print e

	FILE = ''
	dirList=os.listdir(PATH)
	for file in dirList:
		if fnmatch.fnmatch(file, '*mhd.'):
			print PATH+'/'+file

	#TEMP2 = '/data/research/Fiji.app/fiji-linux64 --headless -macro '+MACRO_PATH+'/fijiMacroNRRDSave.ijm '+FILE_GFP+' -batch'
	#TEMP3 = subprocess.Popen(TEMP2, shell=True)
	#print 'BackSubstraction of '+TEMP2
	#TEMP3.communicate()