#!/usr/bin/env python
# -*- coding: utf-8 -*-

import shutil
import fnmatch
import os

# ---------------------------------------------------------------------------------------------------------------------------------------
# Create Folder
# ---------------------------------------------------------------------------------------------------------------------------------------

def main( LOCAL_DATASET_PATH_DATA, PATTERN ):

	FILE = ''
	ENCONTRADO = 0
	dirListGlo=os.listdir(LOCAL_DATASET_PATH_DATA)
	for fileGlo in dirListGlo:
		if fnmatch.fnmatch(fileGlo, PATTERN+'nrrd'):
			FILE = LOCAL_DATASET_PATH_DATA+'/'+fileGlo
			FILE = FILE.rsplit('.',1)[0]
			ENCONTRADO = 1
	if( ENCONTRADO == 1 ):
		return FILE
	else:
		print "ERROR NO SE ENCONTRO EL ARCHIVO"
				#print 'copying: '+GLOBAL_DATASET_PATH+'/'+fileGlo+', to: '+LOCAL_DATASET_PATH_DATA+'/'+fileGlo
				##print(fileGlo)
				#shutil.copy(GLOBAL_DATASET_PATH+'/'+fileGlo, LOCAL_DATASET_PATH_DATA+'/'+fileGlo)
				#FILE = LOCAL_DATASET_PATH_DATA+'/'+fileGlo
				#Primero = 0
			
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