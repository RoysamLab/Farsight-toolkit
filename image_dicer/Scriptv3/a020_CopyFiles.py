#!/usr/bin/env python
# -*- coding: utf-8 -*-

import shutil
import fnmatch
import os
import subprocess

# ---------------------------------------------------------------------------------------------------------------------------------------
# Create Folder
# ---------------------------------------------------------------------------------------------------------------------------------------

def main( GLOBAL_DATASET_PATH, LOCAL_DATASET_PATH_DATA, PATTERN ):

	p = os.uname()
	Primero = 1
	FILE = ''
	dirListGlo=os.listdir(GLOBAL_DATASET_PATH)
	dirListGlo = sorted(dirListGlo)
	for fileGlo in dirListGlo:
		print fileGlo
		if( Primero != 0):
			if fnmatch.fnmatch(fileGlo, PATTERN+'nrrd'):
				dirListlocal=os.listdir(LOCAL_DATASET_PATH_DATA)
				exist = 0
				for fileLoc in dirListlocal:
					if fnmatch.fnmatch(fileLoc, fileGlo):
						exist = 1
				if( exist == 1):
					FILE = LOCAL_DATASET_PATH_DATA+'/'+fileGlo
					FILE = FILE.rsplit('.',1)[0]
					return FILE
				else:
					if( p[1] == 'Farsight-05.EE.UH.EDU' ):
						TEMP2 = 'scp far-05:'+GLOBAL_DATASET_PATH+'/'+fileGlo+' '+LOCAL_DATASET_PATH_DATA+'/'+fileGlo
					if( p[1] == 'Farsight-04.EE.UH.EDU' ):
						TEMP2 = 'scp far-04:'+GLOBAL_DATASET_PATH+'/'+fileGlo+' '+LOCAL_DATASET_PATH_DATA+'/'+fileGlo
					print '\t'+TEMP2
					TEMP8 = subprocess.Popen(TEMP2, shell=True)
					TEMP8.communicate()

					#print 'copying: '+GLOBAL_DATASET_PATH+'/'+fileGlo+', to: '+LOCAL_DATASET_PATH_DATA+'/'+fileGlo
					##print(fileGlo)
					#shutil.copy(GLOBAL_DATASET_PATH+'/'+fileGlo, LOCAL_DATASET_PATH_DATA+'/'+fileGlo)
					FILE = LOCAL_DATASET_PATH_DATA+'/'+fileGlo
					Primero = 0
			
		#if fnmatch.fnmatch(fileGlo, PATTERN+'raw'):
			#print 'copying: '+GLOBAL_DATASET_PATH+'/'+fileGlo+', to: '+LOCAL_DATASET_PATH_DATA+'/'+fileGlo
			##print(fileGlo)
			#shutil.copy(GLOBAL_DATASET_PATH+'/'+fileGlo, LOCAL_DATASET_PATH_DATA+'/'+fileGlo)
		#if fnmatch.fnmatch(fileGlo, '*TRITCdsu.mhd'):
			#print 'copying: '+GLOBAL_DATASET_PATH+'/'+fileGlo+', to: '+LOCAL_DATASET_PATH_DATA+'/'+fileGlo
			##print(fileGlo)
			#shutil.copy(GLOBAL_DATASET_PATH+'/'+fileGlo, LOCAL_DATASET_PATH_DATA+'/'+fileGlo)
		#if fnmatch.fnmatch(fileGlo, '*TRITCdsu.raw'):
			#print 'copying: '+GLOBAL_DATASET_PATH+'/'+fileGlo+', to: '+LOCAL_DATASET_PATH_DATA+'/'+fileGlo
			##print(fileGlo)
			#shutil.copy(GLOBAL_DATASET_PATH+'/'+fileGlo, LOCAL_DATASET_PATH_DATA+'/'+fileGlo)
		#if fnmatch.fnmatch(fileGlo, '*GFPdsu.mhd'):
			#print 'copying: '+GLOBAL_DATASET_PATH+'/'+fileGlo+', to: '+LOCAL_DATASET_PATH_DATA+'/'+fileGlo
			##print(fileGlo)
			#shutil.copy(GLOBAL_DATASET_PATH+'/'+fileGlo, LOCAL_DATASET_PATH_DATA+'/'+fileGlo)
		#if fnmatch.fnmatch(fileGlo, '*GFPdsu.raw'):
			#print 'copying: '+GLOBAL_DATASET_PATH+'/'+fileGlo+', to: '+LOCAL_DATASET_PATH_DATA+'/'+fileGlo
			##print(fileGlo)
			#shutil.copy(GLOBAL_DATASET_PATH+'/'+fileGlo, LOCAL_DATASET_PATH_DATA+'/'+fileGlo)
		#if fnmatch.fnmatch(fileGlo, '*DAPIdsu.mhd'):
			#print 'copying: '+GLOBAL_DATASET_PATH+'/'+fileGlo+', to: '+LOCAL_DATASET_PATH_DATA+'/'+fileGlo
			##print(fileGlo)
			#shutil.copy(GLOBAL_DATASET_PATH+'/'+fileGlo, LOCAL_DATASET_PATH_DATA+'/'+fileGlo)
		#if fnmatch.fnmatch(fileGlo, '*DAPIdsu.raw'):
			#print 'copying: '+GLOBAL_DATASET_PATH+'/'+fileGlo+', to: '+LOCAL_DATASET_PATH_DATA+'/'+fileGlo
			##print(fileGlo)
			#shutil.copy(GLOBAL_DATASET_PATH+'/'+fileGlo, LOCAL_DATASET_PATH_DATA+'/'+fileGlo)
	if( FILE == '' ):
		return FILE
	FILE = FILE.rsplit('.',1)[0]
	return FILE
	

if __name__ == "__main__":
	main()