#!/usr/bin/env python
# -*- coding: utf-8 -*-

import shutil
import fnmatch
import os
import subprocess
import os.path

# ---------------------------------------------------------------------------------------------------------------------------------------
# Create Folder
# ---------------------------------------------------------------------------------------------------------------------------------------

def main( FARSIGHT_BIN_EXE, FILE_Cy5_MED, FILE_Cy5, IMAGE_TYPE ):

	#if os.path.exists(FILE_Cy5_MED+'.nrrd'):
	#	print "Median filter already exist"
	#else:
	print "Median Filter does not exist"
	#runCopy_db_log = LOCAL_DATASET_PATH_LOG +'/runCopyProjections.log'
	#TEMP = FARSIGHT_BIN_EXE+'/ftkMainDarpa PROJECTION '+FILE_GFP+' '+LOCAL_DATASET_PATH_DATA_DEBUG+' > '+runCopy_db_log+' 2>&1'
	print FILE_Cy5
	print FILE_Cy5_MED
	TEMP = FARSIGHT_BIN_EXE+'/ftkMainDarpa MEDIAN '+FILE_Cy5+'.nrrd '+FILE_Cy5_MED+'.nrrd '+IMAGE_TYPE
	TEMP2 = subprocess.Popen(TEMP, shell=True)
	print 'Median Filtering of '+FILE_Cy5
	TEMP2.communicate()
	#TEMP_FILE = open(runDistMap_log, 'a')
	#TEMP_FILE.write('\nCOMMAND: '+TEMP+'\n')
	#TEMP_FILE.close()

if __name__ == "__main__":
	main( FARSIGHT_BIN_EXE, FILE_Cy5_MED, FILE_Cy5, IMAGE_TYPE )
