#!/usr/bin/env python
# -*- coding: utf-8 -*-

import shutil
import fnmatch
import os

# ---------------------------------------------------------------------------------------------------------------------------------------
# Create Folder
# ---------------------------------------------------------------------------------------------------------------------------------------

def main( LOCAL_DEB_DATASET_PATH, LOCAL_DATASET_PATH_EXE, LOCAL_DATASET_PATH_LOG, LOCAL_DATASET_PATH_DEBUG, LOCAL_DATASET_PATH_DATA, LOCAL_DATASET_PATH_DATA_DEBUG, LOCAL_DATASET_PATH_SEGM, LOCAL_DATASET_PATH_SEGM_DEBUG, LOCAL_DATASET_PATH_SEGM_DEBUG_L2, LOCAL_DATASET_PATH_SEGM_TEMP, LOCAL_DATASET_PATH_TRAC, LOCAL_DATASET_PATH_TRAC_DEBUG, LOCAL_DATASET_PATH_TRAC_DEBUG_L2, LOCAL_DATASET_PATH_TRAC_TEMP, GLOBAL_DATASET_PATH_RESULTS, LOCAL_DATASET_PATH_TRAC_RESULTS, LOCAL_DATASET_PATH_ASTRO_TRAC, LOCAL_DATASET_PATH_ASTRO_TRAC_RESULTS, LOCAL_DATASET_PATH_ASTRO_TRAC_DEBUG, LOCAL_DATASET_PATH_ASTRO_TRAC_DEBUG_L2, LOCAL_DATASET_PATH_ASTRO_TRAC_TEMP ):

	if not os.path.isdir(LOCAL_DEB_DATASET_PATH):
		print 'creating folder: '+LOCAL_DEB_DATASET_PATH
		os.makedirs(LOCAL_DEB_DATASET_PATH)
	else:
		print 'erasing folder: '+LOCAL_DEB_DATASET_PATH
		#for the_file in os.listdir(LOCAL_DEB_DATASET_PATH):
			#file_path = os.path.join(LOCAL_DEB_DATASET_PATH, the_file)
			#try:
				#os.unlink(file_path)
			#except Exception, e:
				#print e

	if not os.path.isdir(LOCAL_DATASET_PATH_EXE):
		print 'creating folder: '+LOCAL_DATASET_PATH_EXE
		os.makedirs(LOCAL_DATASET_PATH_EXE)
	else:
		print 'erasing folder: '+LOCAL_DATASET_PATH_EXE
		for the_file in os.listdir(LOCAL_DATASET_PATH_EXE):
			file_path = os.path.join(LOCAL_DATASET_PATH_EXE, the_file)
			try:
				os.unlink(file_path)
			except Exception, e:
				print e

	if not os.path.isdir(LOCAL_DATASET_PATH_LOG):
		print 'creating folder: '+LOCAL_DATASET_PATH_LOG
		os.makedirs(LOCAL_DATASET_PATH_LOG)
	else:
		print 'erasing folder: '+LOCAL_DATASET_PATH_LOG
		#for the_file in os.listdir(LOCAL_DATASET_PATH_LOG):
			#file_path = os.path.join(LOCAL_DATASET_PATH_LOG, the_file)
			#try:
				#os.unlink(file_path)
			#except Exception, e:
				#print e

	if not os.path.isdir(LOCAL_DATASET_PATH_DEBUG):
		print 'creating folder: '+LOCAL_DATASET_PATH_DEBUG
		os.makedirs(LOCAL_DATASET_PATH_DEBUG)
	else:
		print 'erasing folder: '+LOCAL_DATASET_PATH_DEBUG
		#for the_file in os.listdir(LOCAL_DATASET_PATH_DEBUG):
			#file_path = os.path.join(LOCAL_DATASET_PATH_DEBUG, the_file)
			#try:
				#os.unlink(file_path)
			#except Exception, e:
				#print e

	if not os.path.isdir(LOCAL_DATASET_PATH_DATA):
		print 'creating folder: '+LOCAL_DATASET_PATH_DATA
		os.makedirs(LOCAL_DATASET_PATH_DATA)
	else:
		print 'erasing folder: '+LOCAL_DATASET_PATH_DATA
		#for the_file in os.listdir(LOCAL_DATASET_PATH_DATA):
			#file_path = os.path.join(LOCAL_DATASET_PATH_DATA, the_file)
			#try:
				#os.unlink(file_path)
			#except Exception, e:
				#print e

	if not os.path.isdir(LOCAL_DATASET_PATH_DATA_DEBUG):
		print 'creating folder: '+LOCAL_DATASET_PATH_DATA_DEBUG
		os.makedirs(LOCAL_DATASET_PATH_DATA_DEBUG)
	else:
		print 'erasing folder: '+LOCAL_DATASET_PATH_DATA_DEBUG
		#for the_file in os.listdir(LOCAL_DATASET_PATH_DATA_DEBUG):
			#file_path = os.path.join(LOCAL_DATASET_PATH_DATA_DEBUG, the_file)
			#try:
				#os.unlink(file_path)
			#except Exception, e:
				#print e
 
	if not os.path.isdir(LOCAL_DATASET_PATH_SEGM):
		print 'creating folder: '+LOCAL_DATASET_PATH_SEGM
		os.makedirs(LOCAL_DATASET_PATH_SEGM)
	else:
		print 'erasing folder: '+LOCAL_DATASET_PATH_SEGM
		#for the_file in os.listdir(LOCAL_DATASET_PATH_SEGM):
			#file_path = os.path.join(LOCAL_DATASET_PATH_SEGM, the_file)
			#try:
				#os.unlink(file_path)
			#except Exception, e:
				#print e

	if not os.path.isdir(LOCAL_DATASET_PATH_SEGM_DEBUG):
		print 'creating folder: '+LOCAL_DATASET_PATH_SEGM_DEBUG
		os.makedirs(LOCAL_DATASET_PATH_SEGM_DEBUG)
	else:
		print 'erasing folder: '+LOCAL_DATASET_PATH_SEGM_DEBUG
		#for the_file in os.listdir(LOCAL_DATASET_PATH_SEGM_DEBUG):
			#file_path = os.path.join(LOCAL_DATASET_PATH_SEGM_DEBUG, the_file)
			#try:
				#os.unlink(file_path)
			#except Exception, e:
				#print e

	if not os.path.isdir(LOCAL_DATASET_PATH_SEGM_DEBUG_L2):
		print 'creating folder: '+LOCAL_DATASET_PATH_SEGM_DEBUG_L2
		os.makedirs(LOCAL_DATASET_PATH_SEGM_DEBUG_L2)
	else:
		print 'erasing folder: '+LOCAL_DATASET_PATH_SEGM_DEBUG_L2
		#for the_file in os.listdir(LOCAL_DATASET_PATH_SEGM_DEBUG_L2):
			#file_path = os.path.join(LOCAL_DATASET_PATH_SEGM_DEBUG_L2, the_file)
			#try:
				#os.unlink(file_path)
			#except Exception, e:
				#print e

	if not os.path.isdir(LOCAL_DATASET_PATH_SEGM_TEMP):
		print 'creating folder: '+LOCAL_DATASET_PATH_SEGM_TEMP
		os.makedirs(LOCAL_DATASET_PATH_SEGM_TEMP)
	else:
		print 'erasing folder: '+LOCAL_DATASET_PATH_SEGM_TEMP
		#for the_file in os.listdir(LOCAL_DATASET_PATH_SEGM_TEMP):
			#file_path = os.path.join(LOCAL_DATASET_PATH_SEGM_TEMP, the_file)
			#try:
				#os.unlink(file_path)
			#except Exception, e:
				#print e

	if not os.path.isdir(LOCAL_DATASET_PATH_TRAC):
		print 'creating folder: '+LOCAL_DATASET_PATH_TRAC
		os.makedirs(LOCAL_DATASET_PATH_TRAC)
	else:
		print 'erasing folder: '+LOCAL_DATASET_PATH_TRAC
		#for the_file in os.listdir(LOCAL_DATASET_PATH_TRAC):
			#file_path = os.path.join(LOCAL_DATASET_PATH_TRAC, the_file)
			#try:
				#os.unlink(file_path)
			#except Exception, e:
				#print e
	if not os.path.isdir(LOCAL_DATASET_PATH_TRAC_DEBUG):
		print 'creating folder: '+LOCAL_DATASET_PATH_TRAC_DEBUG
		os.makedirs(LOCAL_DATASET_PATH_TRAC_DEBUG)
	else:
		print 'erasing folder: '+LOCAL_DATASET_PATH_TRAC_DEBUG
		#for the_file in os.listdir(LOCAL_DATASET_PATH_TRAC_DEBUG):
			#file_path = os.path.join(LOCAL_DATASET_PATH_TRAC_DEBUG, the_file)
			#try:
				#os.unlink(file_path)
			#except Exception, e:
				#print e

	if not os.path.isdir(LOCAL_DATASET_PATH_TRAC_DEBUG_L2):
		print 'creating folder: '+LOCAL_DATASET_PATH_TRAC_DEBUG_L2
		os.makedirs(LOCAL_DATASET_PATH_TRAC_DEBUG_L2)
	else:
		print 'erasing folder: '+LOCAL_DATASET_PATH_TRAC_DEBUG_L2
		#for the_file in os.listdir(LOCAL_DATASET_PATH_TRAC_DEBUG_L2):
			#file_path = os.path.join(LOCAL_DATASET_PATH_TRAC_DEBUG_L2, the_file)
			#try:
				#os.unlink(file_path)
			#except Exception, e:
				#print e

	if not os.path.isdir(LOCAL_DATASET_PATH_TRAC_TEMP):
		print 'creating folder: '+LOCAL_DATASET_PATH_TRAC_TEMP
		os.makedirs(LOCAL_DATASET_PATH_TRAC_TEMP)
	else:
		print 'erasing folder: '+LOCAL_DATASET_PATH_TRAC_TEMP
		#for the_file in os.listdir(LOCAL_DATASET_PATH_TRAC_TEMP):
			#file_path = os.path.join(LOCAL_DATASET_PATH_TRAC_TEMP, the_file)
			#try:
				#os.unlink(file_path)
			#except Exception, e:
				#print e

	if not os.path.isdir(GLOBAL_DATASET_PATH_RESULTS):
		print 'creating folder: '+GLOBAL_DATASET_PATH_RESULTS
		os.makedirs(GLOBAL_DATASET_PATH_RESULTS)
	else:
		print 'erasing folder: '+GLOBAL_DATASET_PATH_RESULTS
		#for the_file in os.listdir(GLOBAL_DATASET_PATH_RESULTS):
			#file_path = os.path.join(GLOBAL_DATASET_PATH_RESULTS, the_file)
			#try:
				#os.unlink(file_path)
			#except Exception, e:
				#print e

	if not os.path.isdir(LOCAL_DATASET_PATH_TRAC_RESULTS):
		print 'creating folder: '+LOCAL_DATASET_PATH_TRAC_RESULTS
		os.makedirs(LOCAL_DATASET_PATH_TRAC_RESULTS)
	else:
		print 'erasing folder: '+LOCAL_DATASET_PATH_TRAC_RESULTS
		#for the_file in os.listdir(LOCAL_DATASET_PATH_TRAC_RESULTS):
			#file_path = os.path.join(LOCAL_DATASET_PATH_TRAC_RESULTS, the_file)
			#try:
				#os.unlink(file_path)
			#except Exception, e:
				#print e
				
	if not os.path.isdir(LOCAL_DATASET_PATH_ASTRO_TRAC):
		print 'creating folder: '+LOCAL_DATASET_PATH_ASTRO_TRAC
		os.makedirs(LOCAL_DATASET_PATH_ASTRO_TRAC)
	else:
		print 'erasing folder: '+LOCAL_DATASET_PATH_ASTRO_TRAC
		#for the_file in os.listdir(LOCAL_DATASET_PATH_ASTRO_TRAC):
			#file_path = os.path.join(LOCAL_DATASET_PATH_ASTRO_TRAC, the_file)
			#try:
				#os.unlink(file_path)
			#except Exception, e:
				#print e
				
	if not os.path.isdir(LOCAL_DATASET_PATH_ASTRO_TRAC_RESULTS):
		print 'creating folder: '+LOCAL_DATASET_PATH_ASTRO_TRAC_RESULTS
		os.makedirs(LOCAL_DATASET_PATH_ASTRO_TRAC_RESULTS)
	else:
		print 'erasing folder: '+LOCAL_DATASET_PATH_ASTRO_TRAC_RESULTS
		#for the_file in os.listdir(LOCAL_DATASET_PATH_ASTRO_TRAC_RESULTS):
			#file_path = os.path.join(LOCAL_DATASET_PATH_ASTRO_TRAC_RESULTS, the_file)
			#try:
				#os.unlink(file_path)
			#except Exception, e:
				#print e
				
	if not os.path.isdir(LOCAL_DATASET_PATH_ASTRO_TRAC_DEBUG):
		print 'creating folder: '+LOCAL_DATASET_PATH_ASTRO_TRAC_DEBUG
		os.makedirs(LOCAL_DATASET_PATH_ASTRO_TRAC_DEBUG)
	else:
		print 'erasing folder: '+LOCAL_DATASET_PATH_ASTRO_TRAC_DEBUG
		#for the_file in os.listdir(LOCAL_DATASET_PATH_ASTRO_TRAC_DEBUG):
			#file_path = os.path.join(LOCAL_DATASET_PATH_ASTRO_TRAC_DEBUG, the_file)
			#try:
				#os.unlink(file_path)
			#except Exception, e:
				#print e
				
	if not os.path.isdir(LOCAL_DATASET_PATH_ASTRO_TRAC_DEBUG_L2):
		print 'creating folder: '+LOCAL_DATASET_PATH_ASTRO_TRAC_DEBUG_L2
		os.makedirs(LOCAL_DATASET_PATH_ASTRO_TRAC_DEBUG_L2)
	else:
		print 'erasing folder: '+LOCAL_DATASET_PATH_ASTRO_TRAC_DEBUG_L2
		#for the_file in os.listdir(LOCAL_DATASET_PATH_ASTRO_TRAC_DEBUG_L2):
			#file_path = os.path.join(LOCAL_DATASET_PATH_ASTRO_TRAC_DEBUG_L2, the_file)
			#try:
				#os.unlink(file_path)
			#except Exception, e:
				#print e
				
	if not os.path.isdir(LOCAL_DATASET_PATH_ASTRO_TRAC_TEMP):
		print 'creating folder: '+LOCAL_DATASET_PATH_ASTRO_TRAC_TEMP
		os.makedirs(LOCAL_DATASET_PATH_ASTRO_TRAC_TEMP)
	else:
		print 'erasing folder: '+LOCAL_DATASET_PATH_ASTRO_TRAC_TEMP
		#for the_file in os.listdir(LOCAL_DATASET_PATH_ASTRO_TRAC_TEMP):
			#file_path = os.path.join(LOCAL_DATASET_PATH_ASTRO_TRAC_TEMP, the_file)
			#try:
				#os.unlink(file_path)
			#except Exception, e:
				#print e

if __name__ == "__main__":
	main()