#!/usr/bin/env python
# -*- coding: utf-8 -*-

import shutil
import fnmatch
import os
import subprocess
import os.path
import time

MAIN_PATH = '/FSdata/data/DARPA_MOSAICS'
#DATA_FOLDER = ['/0131','/0113','/0117','/0120','/0123','/0128','/0323','/0405'] # With fiji
DATA_FOLDER = ['/0409','/0410','/0412','/1206'] # With c++
#DATA_FOLDER = ['/0117_test']

#MAIN_PATH = '/FSdata/data'
#DATA_FOLDER = ['/0117_Tile_mhd_smaller','/0117_Tile_mhd']

LOCAL_PATH = '/data/nicolas/dataNew/TEMP'

FARSIGHT_BIN_EXE = '/data/nicolas/farsight_updated/bin/exe'

print "# Create Folder"

for fileName in DATA_FOLDER:
	PATH = MAIN_PATH+fileName
	print "-------------------------------------------------------------------------------------------------"
	print " FOLDER: "+fileName
	print "-------------------------------------------------------------------------------------------------"
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

	#print '\t'+PATH
	#ps = []
	startCopy = time.clock()
	fileNameInDirMhd = ''
	fileNameInDirRaw = ''
	dirList=os.listdir(PATH)
	for fileNameInDir in dirList:
		if fnmatch.fnmatch(fileNameInDir, '*.mhd'):
			fileNameInDirMhd = LOCAL_PATH+'/'+fileNameInDir
			shutil.copy(PATH+'/'+fileNameInDir, fileNameInDirMhd)
			fileNameInDirRaw = LOCAL_PATH+'/'+fileNameInDir.rsplit('.',1)[0]+'.raw'
			fileNameInDirRawGlob = PATH+'/'+fileNameInDir.rsplit('.',1)[0]+'.raw'

			shutil.copy(fileNameInDirRawGlob,fileNameInDirRaw)
			
		#if fnmatch.fnmatch(fileNameInDir, '*.raw'):
			#fileNameInDirRaw = LOCAL_PATH+'/'+fileNameInDir
			#shutil.copy(PATH+'/'+fileNameInDir, fileNameInDirRaw)

	#elapsedCopy = (time.clock() - startCopy)
	#print "\t\tCOPYALL TIME: "+str(elapsedCopy)
	#print "-----------------------------------"

	##ps = []
	#dirList=os.listdir(LOCAL_PATH)
	#for fileNameInDir in dirList:
		#if fnmatch.fnmatch(fileNameInDir, '*.mhd'):
			#startSave = time.clock()
			FILE = LOCAL_PATH+'/'+fileNameInDir.rsplit('.',1)[0]
			#TEMP2 = '/data/research/Fiji.app/fiji-linux64 --headless -macro '+LOCAL_PATH+'/fijiMacroNRRD.ijm '+FILE+' -batch'
			TEMP2 = FARSIGHT_BIN_EXE+'/ftkMainDarpa SAVENRRD '+FILE+'.mhd '
			TEMP3 = subprocess.Popen(TEMP2, shell=True)
			#ps.append(TEMP3)
			#print 'NRRD of '+TEMP2
			TEMP3.communicate()
			##elapsedSave = (time.clock() - startSave)
			##print "\t\tSAVE TIME: "+str(elapsedSave)

	#for p in ps:
		#p.wait()

	#dirList=os.listdir(LOCAL_PATH)
	#for fileNameInDir in dirList:
		#if fnmatch.fnmatch(fileNameInDir, '*.mhd'):
			startCopy = time.clock()
			print FILE+'.nrrd --->>> to --->>> '+PATH+'_NRRD'
			shutil.copy(FILE+'.nrrd', PATH+'_NRRD')
			elapsedCopy = (time.clock() - startCopy)
			print "\t\tCOPY TIME: "+str(elapsedCopy)

			os.remove(FILE+'.mhd')
			os.remove(FILE+'.raw')
			os.remove(FILE+'.nrrd')
