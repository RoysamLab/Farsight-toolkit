#!/usr/bin/env python
# -*- coding: utf-8 -*-

import shutil
import fnmatch
import os
import subprocess
import os.path
import shutil
import filecmp
#import md5
import hashlib
import datetime

# ---------------------------------------------------------------------------------------------------------------------------------------
# Create Folder
# ---------------------------------------------------------------------------------------------------------------------------------------

def main( LOCAL_DATASET_PATH_DATA, REMOVE_MONTAGES, GLOBAL_DATASET_PATH_RESULTS ):

	FILE = ''
	dirListGlo=os.listdir(LOCAL_DATASET_PATH_DATA)
	TEMP_FILE = open(LOCAL_DATASET_PATH_DATA+'/checksums', 'w')
	TEMP_FILE.write('CheckSums\n')
	dirListGlo = sorted(dirListGlo)
	for fileGlo in dirListGlo:
		print '\t'+str(datetime.datetime.now())
		FILE = LOCAL_DATASET_PATH_DATA+'/'+fileGlo
		f1 = file(FILE,'rb')
		#SUM = md5.new(f1.read()).hexdigest()
		md5 = hashlib.md5()
		md5.update(f1.read())
		SUM = md5.hexdigest()
		TEMP_FILE.write(SUM+'\t'+fileGlo+'\n')
		print "Check sum of: "+fileGlo+", done. Now It will be send ",
		if( REMOVE_MONTAGES == 1 ):
			TEMP2 = 'scp -C '+FILE+' far-01:'+GLOBAL_DATASET_PATH_RESULTS+'/'
			print '\t'
			TEMP2 = TEMP2 #+';echo "------------------->>>>>>>>>>>>>>>>>>> FINCOPYINGDATA <<<<<<<<<<<<<<<<<<<<-----------------"OF: '+DATA_FOLDER+", START: "+STARTOFCOPY
			TEMP3 = subprocess.Popen(TEMP2, shell=True)
			TEMP3.communicate()
			TEMP2 = 'rm '+FILE
			TEMP3 = subprocess.Popen(TEMP2, shell=True)
			TEMP3.communicate()
			print ' Send and removal done.'
		else:
			TEMP2 = 'scp -C '+FILE+' far-01:'+GLOBAL_DATASET_PATH_RESULTS+'/'
			print '\t'+TEMP2
			TEMP3 = subprocess.Popen(TEMP2, shell=True)
			TEMP3.communicate()
			print ' Send done.'

	TEMP_FILE.close()


if __name__ == "__main__":
	main(  )