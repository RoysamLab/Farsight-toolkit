#!/usr/bin/env python
# -*- coding: utf-8 -*-

import shutil
import fnmatch
import os
import subprocess

# ---------------------------------------------------------------------------------------------------------------------------------------
# Create Folder
# ---------------------------------------------------------------------------------------------------------------------------------------

def main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_DATA_DEBUG, FILE_GFP, runCopy_db_log, PROJEOPT, IMAGETYPE ):
	
	NEW_NAME = LOCAL_DATASET_PATH_DATA_DEBUG+'/'+os.path.basename(FILE_GFP)
	print "\t\t"+NEW_NAME
	if( os.path.exists(NEW_NAME+'zPro_X.tif') & os.path.exists(NEW_NAME+'zPro_Y.tif') & os.path.exists(NEW_NAME+'zPro_Z.tif') & os.path.exists(NEW_NAME+'zPro_X_Re.tif') & os.path.exists(NEW_NAME+'zPro_Y_Re.tif') & os.path.exists(NEW_NAME+'zPro_Z_Re.tif') ):
		print "Projection already exist"
	else:
		print "Projection does not exist"
		#runCopy_db_log = LOCAL_DATASET_PATH_LOG +'/runCopyProjections.log'
		TEMP = FARSIGHT_BIN_EXE+'/ftkMainDarpa PROJECTION '+FILE_GFP+'.nrrd '+LOCAL_DATASET_PATH_DATA_DEBUG+' '+PROJEOPT+' '+IMAGETYPE+'>> '+runCopy_db_log+' 2>&1'
		TEMP2 = subprocess.Popen(TEMP, shell=True)
		print 'Projection of '+FILE_GFP
		TEMP2.communicate()
		TEMP_FILE = open(runCopy_db_log, 'a')
		TEMP_FILE.write('\nCOMMAND: '+TEMP+'\n')
		TEMP_FILE.close()

#FILE = FILE.rsplit('/',1)[0]

if __name__ == "__main__":
	main()