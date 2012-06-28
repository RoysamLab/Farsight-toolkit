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

def main( FARSIGHT_BIN_EXE, optionsTracing, runTrac_log, LOCAL_DATASET_PATH_DATA ):

	# THIS FILE COMPARISON IS NOT WORKIN, IS GIVEIN AND EOF DIFFERENCE
	flag = 0
	sameFile = 0
	print "HERE: "+LOCAL_DATASET_PATH_DATA+'/options_tracing'
	if( os.path.exists(LOCAL_DATASET_PATH_DATA+'/options_tracing') ):
		flag = 1
		print "file exists" + str(filecmp.cmp(optionsTracing, LOCAL_DATASET_PATH_DATA+'/options_tracing'))
	if( (flag == 1)):
		print optionsTracing+" "+LOCAL_DATASET_PATH_DATA+'/options_tracing'
		if( filecmp.cmp(optionsTracing, LOCAL_DATASET_PATH_DATA+'/options_tracing') ):
			sameFile = 1;
			print "file is same"
	if( sameFile == 1 ):
		print 'Tracing result already exists, same parameters'
	else:
		TEMP = FARSIGHT_BIN_EXE+'/ftkMainDarpa TRACE '+optionsTracing+' >> '+runTrac_log+' 2>&1'
		TEMP2 = subprocess.Popen(TEMP, shell=True)
		print 'Tracing'
		TEMP2.communicate()
		TEMP_FILE = open(runTrac_log, 'a')
		TEMP_FILE.write('\nCOMMAND: '+TEMP+'\n')
		TEMP_FILE.close()
		shutil.copy(optionsTracing, LOCAL_DATASET_PATH_DATA+'/options_tracing')


if __name__ == "__main__":
	main( FARSIGHT_BIN_EXE, optionsSegm, runSegm_log )