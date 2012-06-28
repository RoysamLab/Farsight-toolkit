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

def main( FILE_GFP_BS_CV_RE_8bit, FARSIGHT_BIN_EXE, optionsSegm, runSegm_log, LOCAL_DATASET_PATH_DATA ):

	flag = 0
	sameFile = 0
	if( os.path.exists(LOCAL_DATASET_PATH_DATA+'/options_segmentation') ):
		flag = 1
	if( (flag == 1)):
		if( filecmp.cmp(optionsSegm, LOCAL_DATASET_PATH_DATA+'/options_segmentation') ):
			sameFile = 1;
	if( sameFile == 1 ):
		print 'Segmentation result already exists, same parameters'
	else:
		TEMP = FARSIGHT_BIN_EXE+'/ftkMainDarpa SEGMENT '+optionsSegm+' >> '+runSegm_log+' 2>&1'
		TEMP2 = subprocess.Popen(TEMP, shell=True)
		print 'Segmentation'
		TEMP2.communicate()
		TEMP_FILE = open(runSegm_log, 'a')
		TEMP_FILE.write('\nCOMMAND: '+TEMP+'\n')
		TEMP_FILE.close()
		shutil.copy(optionsSegm, LOCAL_DATASET_PATH_DATA+'/options_segmentation')


if __name__ == "__main__":
	main( FARSIGHT_BIN_EXE, optionsSegm, runSegm_log )