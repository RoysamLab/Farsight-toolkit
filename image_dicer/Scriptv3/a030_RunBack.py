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

def main( LOCAL_DATASET_PATH_PARAMETERS, FILE_GFP_BS, FILE_GFP ):
	if os.path.exists(FILE_GFP_BS+'.nrrd'):
		print "GFP BS already exist"
	else:
		print "GFP BS does not exist"
		#runCopy_db_log = LOCAL_DATASET_PATH_LOG +'/runCopyProjections.log'
		#TEMP = FARSIGHT_BIN_EXE+'/ftkMainDarpa PROJECTION '+FILE_GFP+' '+LOCAL_DATASET_PATH_DATA_DEBUG+' > '+runCopy_db_log+' 2>&1'
		TEMP2 = '/data/research/Fiji.app/fiji-linux64 --headless -macro '+LOCAL_DATASET_PATH_PARAMETERS+'/fijiMacroBackSubst.ijm '+FILE_GFP+' -batch'
		TEMP3 = subprocess.Popen(TEMP2, shell=True)
		print 'BackSubstraction of '+TEMP2
		TEMP3.communicate()
	#TEMP_FILE = open(runCopy_db_log, 'a')
	#TEMP_FILE.write('\nCOMMAND: '+TEMP+'\n')
	#TEMP_FILE.close()

if __name__ == "__main__":
	main( LOCAL_DATASET_PATH_PARAMETERS, FILE_GFP_BS, FILE_GFP )