#!/usr/bin/env python
# -*- coding: utf-8 -*-

import shutil
import fnmatch
import os
import subprocess

# ---------------------------------------------------------------------------------------------------------------------------------------
# Create Folder
# ---------------------------------------------------------------------------------------------------------------------------------------

def main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_SEGM_DEBUG, FILE_GFP, FILE_SOMA, NAME, runCopy_db_log ):

	NEW_NAME = LOCAL_DATASET_PATH_SEGM_DEBUG+'/'
	if( os.path.exists(NEW_NAME+'zRGB_'+NAME+'Pro_Z.tif') & os.path.exists(NEW_NAME+'zRGB_'+NAME+'Pro_Y.tif') & os.path.exists(NEW_NAME+'zRGB_'+NAME+'Pro_X.tif') ):
		print "Projection already exist"
	else:
		print "Projection does not exist"
		#runCopy_db_log = LOCAL_DATASET_PATH_LOG +'/runCopyProjections.log'
		TEMP = FARSIGHT_BIN_EXE+'/ftkMainDarpa PROJECTIONRGB '+FILE_GFP+'.nrrd '+FILE_SOMA+'.nrrd '+LOCAL_DATASET_PATH_SEGM_DEBUG+' '+NAME+' >> '+runCopy_db_log+' 2>&1'
		TEMP2 = subprocess.Popen(TEMP, shell=True)
		print 'Projection RGB of '+FILE_GFP
		TEMP2.communicate()
		TEMP_FILE = open(runCopy_db_log, 'a')
		TEMP_FILE.write('\nCOMMAND: '+TEMP+'\n')
		TEMP_FILE.close()

if __name__ == "__main__":
	main()



#a023_ProjectRGB.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_SEGM_DEBUG, FILE_GFP_BS_CV_RE_bit, FILE_GFP_BS_CV_RE_bit+'_soma', runSegm_db_log )