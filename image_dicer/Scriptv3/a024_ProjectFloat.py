#!/usr/bin/env python
# -*- coding: utf-8 -*-

import shutil
import fnmatch
import os
import subprocess

# ---------------------------------------------------------------------------------------------------------------------------------------
# Create Folder
# ---------------------------------------------------------------------------------------------------------------------------------------

def main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_TRAC_DEBUG, FILE_GFP, runCopy_db_log, PROJEOPT, IMAGETYPE ):

	NEW_NAME = LOCAL_DATASET_PATH_TRAC_DEBUG+'/'+os.path.basename(FILE_GFP)
	if( os.path.exists(NEW_NAME+'zPro_Z.nrrd') & os.path.exists(NEW_NAME+'zPro_Y.nrrd') & os.path.exists(NEW_NAME+'zPro_X.nrrd') ):
		print "Projection already exist"
	else:
		print "Projection does not exist"
		#runCopy_db_log = LOCAL_DATASET_PATH_LOG +'/runCopyProjections.log'
		TEMP = FARSIGHT_BIN_EXE+'/ftkMainDarpa PROJECTIONFLO '+FILE_GFP+'.nrrd '+LOCAL_DATASET_PATH_TRAC_DEBUG+' '+PROJEOPT+' '+IMAGETYPE+' >> '+runCopy_db_log+' 2>&1'
		TEMP2 = subprocess.Popen(TEMP, shell=True)
		print 'Projection FLO of '+FILE_GFP
		TEMP2.communicate()
		TEMP_FILE = open(runCopy_db_log, 'a')
		TEMP_FILE.write('\nCOMMAND: '+TEMP+'\n')
		TEMP_FILE.close()

if __name__ == "__main__":
	main()



#a023_ProjectRGB.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_SEGM_DEBUG, FILE_GFP_BS_CV_RE_bit, FILE_GFP_BS_CV_RE_bit+'_soma', runSegm_db_log )