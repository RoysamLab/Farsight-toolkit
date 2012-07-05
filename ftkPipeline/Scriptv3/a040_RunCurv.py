#!/usr/bin/env python
# -*- coding: utf-8 -*-

import shutil
import fnmatch
import os
import subprocess
import os.path
import shutil
import filecmp

# ---------------------------------------------------------------------------------------------------------------------------------------
# Create Folder
# ---------------------------------------------------------------------------------------------------------------------------------------

def main( FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_PARAMETERS, FILE_GFP_BS, runCurv_db_log, LOCAL_DATASET_PATH_DATA ):

	flag = 0
	sameFile = 0
	if( os.path.exists(LOCAL_DATASET_PATH_DATA+'/options_curvelets') ):
		flag = 1
	if( (flag == 1)):
		if( filecmp.cmp(LOCAL_DATASET_PATH_PARAMETERS+'/options_curvelets', LOCAL_DATASET_PATH_DATA+'/options_curvelets') ):
			sameFile = 1;
	if( sameFile == 1 ):
		print 'Curvelets result already exists, same parameters'
	else:
		TEMP = FARSIGHT_BIN_EXE+'/curvelets '+FILE_GFP_BS+'.nrrd '+LOCAL_DATASET_PATH_PARAMETERS+'/options_curvelets >> '+runCurv_db_log+' 2>&1'
		TEMP3 = subprocess.Popen(TEMP, shell=True)
		print 'Curvelets of '+FILE_GFP_BS
		TEMP3.communicate()
		TEMP_FILE = open(runCurv_db_log, 'a')
		TEMP_FILE.write('\nCOMMAND: '+TEMP+'\n')
		TEMP_FILE.close()
		shutil.copy(LOCAL_DATASET_PATH_PARAMETERS+'/options_curvelets', LOCAL_DATASET_PATH_DATA+'/options_curvelets')

if __name__ == "__main__":
	main( FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_PARAMETERS, FILE_GFP_BS, runCurv_db_log )