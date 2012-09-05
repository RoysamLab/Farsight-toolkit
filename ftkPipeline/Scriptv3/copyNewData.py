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

def main(  ):


	REMOTE_PATH = "/media/BIA\ LAB\ 03_/Susie"
	LOCAL_PATH = "/FSdata/data/DARPA_RAW/Susie"

	if not os.path.isdir(LOCAL_PATH):
		print 'creating folder: '+LOCAL_PATH
		os.makedirs(LOCAL_PATH)
	else:
		print 'erasing folder: '+LOCAL_PATH
		#for the_file in os.listdir(LOCAL_DEB_DATASET_PATH):
			#file_path = os.path.join(LOCAL_DEB_DATASET_PATH, the_file)
			#try:
				#os.unlink(file_path)
			#except Exception, e:
				#print e

	dirListRemo=os.listdir(REMOTE_PATH)
	dirListRemo = sorted(dirListRemo)
	for fileRemo in dirListRemo:
		print '\t'+str(datetime.datetime.now())
		FILE_REMO = REMOTE_PATH+'/'+fileRemo
		# Checksum Original
		f1 = file(FILE_REMO,'rb')
		md5 = hashlib.md5()
		md5.update(f1.read())
		SUM_ORIGINAL = md5.hexdigest()
		print "Check sum of: "+fileRemo+", done. Now It will be send ",
		# Send
		TEMP2 = 'scp -C '+FILE_REMO+' far-01:'+LOCAL_PATH
		print TEMP2
		TEMP3 = subprocess.Popen(TEMP2, shell=True)
		TEMP3.communicate();
		# Checksum after copy
		FILE_LOCAL = LOCAL_PATH+'/'+fileRemo
		f1 = file(FILE_LOCAL,'rb')
		md5 = hashlib.md5()
		md5.update(f1.read())
		SUM_NEW = md5.hexdigest()
		# Compare the sums
		if( SUM_ORIGINAL != SUM_NEW )
			print "\t\tERROR THE CHECK SUM OF "+fileRemo+" IS NOT THE SAME"
			


if __name__ == "__main__":
	main(  )
