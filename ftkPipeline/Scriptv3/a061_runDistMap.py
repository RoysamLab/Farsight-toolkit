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

def main( FARSIGHT_BIN_EXE, FILE_LABEL_DM, FILE_LABEL, runDistMap_log ):

	if os.path.exists(FILE_LABEL_DM+'.nrrd'):
		print "Distance map already exist"
	else:
		print "Distance map does not exist"
		#runCopy_db_log = LOCAL_DATASET_PATH_LOG +'/runCopyProjections.log'
		#TEMP = FARSIGHT_BIN_EXE+'/ftkMainDarpa PROJECTION '+FILE_GFP+' '+LOCAL_DATASET_PATH_DATA_DEBUG+' > '+runCopy_db_log+' 2>&1'
		TEMP = FARSIGHT_BIN_EXE+'/ftkMainDarpa DISTANCE_MAP '+FILE_LABEL+'.nrrd '+FILE_LABEL_DM+'.nrrd'+' >> '+runDistMap_log+' 2>&1'
		TEMP2 = subprocess.Popen(TEMP, shell=True)
		print 'Distance Map of '+FILE_LABEL
		TEMP2.communicate()
		TEMP_FILE = open(runDistMap_log, 'a')
		TEMP_FILE.write('\nCOMMAND: '+TEMP+'\n')
		TEMP_FILE.close()

if __name__ == "__main__":
	main( FARSIGHT_BIN_EXE, FILE_LABEL_DM, FILE_LABEL, runDistMap_log )