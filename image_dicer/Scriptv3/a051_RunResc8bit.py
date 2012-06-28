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

def main( FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_PARAMETERS, OUTPUT, INPUT, runRescale_log ):

	if os.path.exists(OUTPUT+'.nrrd'):
		print "Rescale Exits already exist"
	else:
		print "Rescale does not exist"
		#runCopy_db_log = LOCAL_DATASET_PATH_LOG +'/runCopyProjections.log'
		#TEMP = FARSIGHT_BIN_EXE+'/ftkMainDarpa PROJECTION '+FILE_GFP+' '+LOCAL_DATASET_PATH_DATA_DEBUG+' > '+runCopy_db_log+' 2>&1'
		TEMP = FARSIGHT_BIN_EXE+'/ftkMainDarpa RESCALE_8BIT '+INPUT+'.nrrd '+OUTPUT+'.nrrd'+' >> '+runRescale_log+' 2>&1'
		TEMP2 = subprocess.Popen(TEMP, shell=True)
		print 'Rescale of '+INPUT
		TEMP2.communicate()
		TEMP_FILE = open(runRescale_log, 'a')
		TEMP_FILE.write('\nCOMMAND: '+TEMP+'\n')
		TEMP_FILE.close()



if __name__ == "__main__":
	main( FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_PARAMETERS, OUTPUT, INPUT, runRescale_log )