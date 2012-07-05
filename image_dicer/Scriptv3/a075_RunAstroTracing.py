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

def main( FILE_TRI_BS_CV_RE_bit, FARSIGHT_BIN_EXE, optionsAstroTrac, runAstroTrac_log, LOCAL_DATASET_PATH_DATA ):

	flag = 0
	sameFile = 0
	if( os.path.exists(LOCAL_DATASET_PATH_DATA+'/options_astrocyte_tracing') ):
		flag = 1
	if( (flag == 1)):
		if( filecmp.cmp(optionsAstroTrac, LOCAL_DATASET_PATH_DATA+'/options_astrocyte_tracing') ):
			sameFile = 1;
	if( sameFile == 1 ):
		print 'Astrocyte_Tracing result already exists, same parameters'
	else:
		TEMP = FARSIGHT_BIN_EXE+'/ftkMainDarpa ASTRO_TRACE '+optionsAstroTrac+' >> '+runAstroTrac_log+' 2>&1'
		TEMP2 = subprocess.Popen(TEMP, shell=True)
		print 'Astrocyte_Tracing'
		TEMP2.communicate()
		TEMP_FILE = open(runAstroTrac_log, 'a')
		TEMP_FILE.write('\nCOMMAND: '+TEMP+'\n')
		TEMP_FILE.close()
		shutil.copy(optionsAstroTrac, LOCAL_DATASET_PATH_DATA+'/options_astrocyte_tracing')


if __name__ == "__main__":
	main( FARSIGHT_BIN_EXE, optionsAstroTrac, runAstroTrac_log )