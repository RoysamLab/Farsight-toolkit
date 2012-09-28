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

def main( FARSIGHT_BIN_EXE, FILE_DEVICE, INPUT_TABLE, OUTPUT_TABLE, IMAGE_TYPE, runDistToDevice_log ):

	if os.path.exists(OUTPUT_TABLE):
		print "Distance to device feature already exist"
	else:
		print "Distance to device feature does not exist"
		#runCopy_db_log = LOCAL_DATASET_PATH_LOG +'/runCopyProjections.log'
		#TEMP = FARSIGHT_BIN_EXE+'/ftkMainDarpa PROJECTION '+FILE_GFP+' '+LOCAL_DATASET_PATH_DATA_DEBUG+' > '+runCopy_db_log+' 2>&1'
		TEMP = FARSIGHT_BIN_EXE+'/ftkMainDarpa DISTANCE_TO_DEVICE '+FILE_DEVICE+'.nrrd '+INPUT_TABLE+' '+OUTPUT_TABLE+' '+IMAGE_TYPE+' >> '+runDistToDevice_log+' 2>&1'
		TEMP2 = subprocess.Popen(TEMP, shell=True)
		print 'Computing distance to device feature for '+INPUT_TABLE
		TEMP2.communicate()
		TEMP_FILE = open(runDistToDevice_log, 'a')
		TEMP_FILE.write('\nCOMMAND: '+TEMP+'\n')
		TEMP_FILE.close()

if __name__ == "__main__":
	main( FARSIGHT_BIN_EXE, FILE_DEVICE, INPUT_TABLE, OUTPUT_TABLE, IMAGE_TYPE, runDistToDevice_log )