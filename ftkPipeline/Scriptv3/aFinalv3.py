#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Check MNT params
# Curv params
# xTile yTile
# Origin and destination folders

# 1 Sync
# 2 Create folder
# 3 setupt parameters
# 4 run

import shutil
import fnmatch
import os
import subprocess
import os.path
import time
import glob
import sys
import inspect
import filecmp
import datetime
import getpass #ask username


import a010_CreateFolders
import a019_FindFile
import a020_CopyFiles
import a021_Project
import a022_Project8bits
import a023_ProjectRGB
import a024_ProjectFloat
import a030_RunBack
import a040_RunCurv
import a050_RunResc
import a051_RunResc8bit
import a052_RunResc8bitMhd
import a060_RunSegmentation
import a070_RunTracing
import a075_RunAstroTracing
#import a080_RunCheckSum
import a080_RunCheckSumNew
#import 
#import a080_ComputeHistogram

#tail -f `ls -tr | tail -n 1`


class MyWriter:
	def __init__(self,stdout, filename):
		self.stdout = stdout
		self.logfile = file(filename, 'w')
	def write(self, text):
		self.stdout.write(text)
		self.logfile.write(text)
		self.logfile.flush()

	def close(self):
		self.stdout.close()
		self.logfile.close()

TEMP = inspect.getfile(inspect.currentframe())
TEMP = os.path.basename(TEMP)
TEMP = os.path.splitext(TEMP)[0]
writerLog = MyWriter(sys.stdout, TEMP+'_log')
sys.stdout = writerLog


TEST_RUN = 0
if len(sys.argv) == 2:
	if(sys.argv[1] =='TEST' ):
		TEST_RUN = 1
	else:
		TEST_RUN = 0


#def main(DATA_FOLDER):

# Parameters
TRY = 07
#DEBUG = 1
#----------------------------------
REMOVE_MONTAGES = 0	# This flag is set in case we want the montages to be removed after the process is done, especially when running many montages in serial we want to make sure not to
MOVE_RESULTS = 0	# If 0 the results will be keep
			# if 1 the results will be moved
			# if 2 the results will be copied (keep and copy to FSDATA)
			# if 3 move everysingle file, exept the folder, which are copied
MOVE_LOCALLY = 0
#----------------------------------
REMOVE_TEMP_SEGM = 1
REMOVE_TEMP_TRAC = 1
#----------------------------------
MOVE_INTER_STEPS = 0

SMALLIMAGE = '1'	# if the image is small

runCopy = 1		# Flag to move images
runCopy_db = 1
runMake = 1		# Flag Make Farsight
runBack = 1		# Flag to run background substraction
runBack_db = 1
runDistMap = 1		# Flag to run distance map computation
runDistMap_db = 1

runCurv = 1		# Flag to run Curvelets
runCurv_db = 1
runRescale = 1
runRescale_db = 1
runRescale_bit = 1
runRescale_bit_db = 1
runSegm = 1		# Flag to run Segmentation
runSegm_db = 1
runTrac = 1		# Flag to run Tracing
runTrac_db = 1
runAstroTrac = 1		# Flag to run Astrocyte Tracing
runAstroTrac_db = 1

#runCurv = 0#1		# Flag to run Curvelets
#runCurv_db = 0#1
#runRescale = 0#1
#runRescale_db = 0#1
#runRescale_bit = 0#1
#runRescale_bit_db = 0#1
#runSegm = 0#1		# Flag to run Segmentation
#runSegm_db = 0#1
#runTrac = 0#1		# Flag to run Tracing
#runTrac_db = 0#1

haveCy5_for = 0
haveTRT_for = 0
haveGFP_for = 1
haveDAP_for = 1

# ---------------------------------------------------------------------------------------------------------------------------------------

SERVER = 'none'
p = os.uname()
if( p[1] == 'Farsight-05.EE.UH.EDU' ):
	SERVER = 'far05'
if( p[1] == 'Farsight-04.EE.UH.EDU' ):
	SERVER = 'far04'
	
	
USERNAME = 'none'
p = getpass.getuser()
if( p[1] == 'vhsomasu' ):
	USERNAME = 'vhsomasu'
if( p[1] == 'nrey' ):
	USERNAME = 'nrey'

#DATA_FOLDER_ALL = ['/0131_test','/0131_test2'] # For testing dont forget the xTile params
#DATA_FOLDER_ALL = ['/0131_test'] # For testing dont forget the xTile params
#DATA_FOLDER_ALL = ['/0113_NRRD_CROPPED']

#DATA_FOLDER_ALL = ['/0113_NRRD_CROPPED','/0117_NRRD_CROPPED','/0120_NRRD_CROPPED','/0123_NRRD_CROPPED','/0128_NRRD_CROPPED','/0131_NRRD_CROPPED','/0323_NRRD_CROPPED','/0405_NRRD_CROPPED','/0409_NRRD_CROPPED','/0410_NRRD_CROPPED','/0412_NRRD_CROPPED','/1206_NRRD_CROPPED']
#DATA_FOLDER_ALL = ['/0113_NRRD_CROPPED','/0117_NRRD_CROPPED','/0120_NRRD_CROPPED','/0123_NRRD_CROPPED','/0128_NRRD_CROPPED','/0131_NRRD_CROPPED']
#DATA_FOLDER_ALL = ['/0323_NRRD_CROPPED','/0405_NRRD_CROPPED','/0409_NRRD_CROPPED','/0410_NRRD_CROPPED','/0412_NRRD_CROPPED','/1206_NRRD_CROPPED']

#DATA_FOLDER_ALL = ['/0113_NRRD','/0117_NRRD','/0120_NRRD','/0123_NRRD','/0128_NRRD','/0131_NRRD','/0323_NRRD','/0405_NRRD','/0409_NRRD','/0410_NRRD','/0412_NRRD','/1206_NRRD']

DATA_FOLDER_ALL = ['/0131_test']
if( SERVER == 'far04' ):
	REMOVE_MONTAGES = 1
	#DATA_FOLDER_ALL = ['/0113_NRRD_CROPPED','/0117_NRRD_CROPPED','/0120_NRRD_CROPPED','/0123_NRRD_CROPPED','/0128_NRRD_CROPPED','/0131_NRRD_CROPPED','/0323_NRRD_CROPPED','/0405_NRRD_CROPPED','/0409_NRRD_CROPPED','/0410_NRRD_CROPPED','/0412_NRRD_CROPPED','/1206_NRRD_CROPPED']
	#DATA_FOLDER_ALL = ['/0113_NRRD','/0117_NRRD','/0120_NRRD','/0123_NRRD','/0128_NRRD','/0131_NRRD','/0323_NRRD','/0405_NRRD','/0409_NRRD','/0410_NRRD','/0412_NRRD','/1206_NRRD']
	#DATA_FOLDER_ALL = ['/0323_NRRD','/0405_NRRD','/0409_NRRD','/0410_NRRD','/0412_NRRD','/1206_NRRD']
	#DATA_FOLDER_ALL = ['/0113_NRRD','/0117_NRRD','/0120_NRRD','/0123_NRRD','/0128_NRRD','/0131_NRRD']
	#DATA_FOLDER_ALL = ['/0117_NRRD','/0120_NRRD','/0123_NRRD','/0128_NRRD','/0131_NRRD']
	#DATA_FOLDER_ALL = ['/0120_NRRD','/0123_NRRD','/0128_NRRD','/0131_NRRD']
	#DATA_FOLDER_ALL = ['/0113_NRRD']
	#DATA_FOLDER_ALL = ['/0131_test']
	#DATA_FOLDER_ALL = ['/0131_test2']
	#DATA_FOLDER_ALL = ['/0131_test3']
	DATA_FOLDER_ALL = ['/0131_test4']
	#DATA_FOLDER_ALL = ['/0120_NRRD']
if( SERVER == 'far05' ):
	REMOVE_MONTAGES = 1
	#DATA_FOLDER_ALL = ['/0113_NRRD','/0117_NRRD','/0120_NRRD','/0123_NRRD','/0128_NRRD','/0131_NRRD','/0323_NRRD','/0405_NRRD','/0409_NRRD','/0410_NRRD','/0412_NRRD','/1206_NRRD']
	DATA_FOLDER_ALL = ['/0323_NRRD','/0405_NRRD','/0409_NRRD','/0410_NRRD','/0412_NRRD','/1206_NRRD']
	#DATA_FOLDER_ALL = ['/0405_NRRD','/0409_NRRD','/0410_NRRD','/0412_NRRD','/1206_NRRD']
	#DATA_FOLDER_ALL = ['/0131_test']
	#DATA_FOLDER_ALL = ['/0131_test2']
	#DATA_FOLDER_ALL = ['/0131_test3']
	#DATA_FOLDER_ALL = ['/0131_test4']
	runMake = 0
if( TEST_RUN==1 ):
	DATA_FOLDER_ALL = ['/0131_test4']
	print "THIS IS A TEST"
	print "THIS IS A TEST"
# ---------------------------------------------------------------------------------------------------------------------------------------

#def main(DATA_FOLDER):

for DATA_FOLDER in DATA_FOLDER_ALL:

	# ---------------------------------------------------------------------------------------------------------------------------------------
	haveCy5 = 1
	haveTRT = 1
	haveGFP = 1
	haveDAP = 1

	if( haveCy5_for == 0 ):
		haveCy5 = 0
	if( haveTRT_for == 0 ):
		haveTRT = 0
	if( haveGFP_for == 0 ):
		haveGFP = 0
	if( haveDAP_for == 0 ):
		haveDAP = 0

	FILE_Cy5 = ''
	FILE_TRI = ''
	FILE_GFP = ''
	FILE_DAP = ''
	FILE_Cy5_BS = ''
	FILE_TRI_BS = ''
	FILE_GFP_BS = ''
	FILE_DAP_BS = ''
	FILE_Cy5_BS_RE = ''
	FILE_TRI_BS_RE = ''
	FILE_GFP_BS_CV_RE = ''
	FILE_DAP_BS_RE = ''
	FILE_Cy5_BS_RE_bit = ''
	FILE_TRI_BS_RE_bit = ''
	FILE_GFP_BS_CV_RE_bit = ''
	FILE_DAP_BS_RE_bit = ''
	
	# ---------------------------------------------------------------------------------------------------------------------------------------
	# Folders names 
	# ---------------------------------------------------------------------------------------------------------------------------------------

	if( USERNAME == 'nrey' ):
		FARSIGHT_BIN = '/data/nicolas/farsight_updated/bin'
		FARSIGHT_BIN_EXE = '/data/nicolas/farsight_updated/bin/exe'

		MAIN_DATA_FOLDER = '/data/nicolas/dataNew'
		MAIN_DEB_DATA_FOLDER = '/data/nicolas/deb'
		
	if( USERNAME = 'vhsomasu' ):
		FARSIGHT_BIN = '/data/vinay/farsight_v1/bin'
		FARSIGHT_BIN_EXE = '/data/vinay/farsight_v1/bin/exe'

		MAIN_DATA_FOLDER = '/data/vinay/dataNew'
		MAIN_DEB_DATA_FOLDER = '/data/vinay/deb'
		

	LOCAL_DATASET_PATH = MAIN_DATA_FOLDER+DATA_FOLDER
	LOCAL_DEB_DATASET_PATH = MAIN_DEB_DATA_FOLDER+'/'+str(TRY)+DATA_FOLDER+'_RESULTS_V3_'+SERVER

	LOCAL_DATASET_PATH_PARAMETERS = LOCAL_DATASET_PATH+'/Parameters'	# ---> This directory has to exist
	LOCAL_DATASET_PATH_EXE = LOCAL_DATASET_PATH+'/Exe'
	LOCAL_DATASET_PATH_LOG = LOCAL_DATASET_PATH+'/Log'
	LOCAL_DATASET_PATH_DEBUG = LOCAL_DATASET_PATH+'/DEBUG'
	#LOCAL_DATASET_PATH_TRACE_SOMASDDIVIDED = LOCAL_DATASET_PATH+'/TracesAndSomasDivided'
	LOCAL_DATASET_PATH_DATA = LOCAL_DATASET_PATH+'/Data'
	LOCAL_DATASET_PATH_DATA_DEBUG = LOCAL_DATASET_PATH_DEBUG+'/Data'
	
	LOCAL_DATASET_PATH_SEGM = LOCAL_DATASET_PATH+'/Segm'
	LOCAL_DATASET_PATH_SEGM_DEBUG = LOCAL_DATASET_PATH_DEBUG+'/Segm'
	LOCAL_DATASET_PATH_SEGM_DEBUG_L2 = LOCAL_DATASET_PATH_SEGM_DEBUG+'/Level2'
	LOCAL_DATASET_PATH_SEGM_TEMP = LOCAL_DATASET_PATH_SEGM+'/Temp'
	
	LOCAL_DATASET_PATH_TRAC = LOCAL_DATASET_PATH+'/Trace'
	LOCAL_DATASET_PATH_TRAC_RESULTS = LOCAL_DATASET_PATH_TRAC+'/Results'
	LOCAL_DATASET_PATH_TRAC_DEBUG = LOCAL_DATASET_PATH_DEBUG+'/Trace'
	LOCAL_DATASET_PATH_TRAC_DEBUG_L2 = LOCAL_DATASET_PATH_TRAC_DEBUG+'/Level2'
	LOCAL_DATASET_PATH_TRAC_TEMP = LOCAL_DATASET_PATH_TRAC+'/Temp'
	
	LOCAL_DATASET_PATH_ASTRO_TRAC = LOCAL_DATASET_PATH+'/Astro_Trace'
	LOCAL_DATASET_PATH_ASTRO_TRAC_RESULTS = LOCAL_DATASET_PATH_ASTRO_TRAC+'/Results'
	LOCAL_DATASET_PATH_ASTRO_TRAC_DEBUG = LOCAL_DATASET_PATH_DEBUG+'/Astro_Trace'
	LOCAL_DATASET_PATH_ASTRO_TRAC_DEBUG_L2 = LOCAL_DATASET_PATH_ASTRO_TRAC_DEBUG+'/Level2'
	LOCAL_DATASET_PATH_ASTRO_TRAC_TEMP = LOCAL_DATASET_PATH_ASTRO_TRAC+'/Temp'

	GLOBAL_DATASET_PATH = "/FSdata/data/DARPA_MOSAICS"+DATA_FOLDER
	GLOBAL_DATASET_PARAMS = "/FSdata/data/DARPA_MOSAICS"+DATA_FOLDER+'_PARAMS_DEVICE_'+'1'#str(TRY)
	#GLOBAL_DATASET_PATH = "/FSdata/data"+DATA_FOLDER
	GLOBAL_DATASET_PATH_RESULTS = '/FSdata/data/DARPA_RESULTS'+'/'+str(TRY)+DATA_FOLDER+'_RESULTS_V3_'+SERVER

	print "# ---------------------------------------------------------------------------------------------------------------------------------------"
	print " Make Farsight: "+DATA_FOLDER+' '+str(datetime.datetime.now())
	print "# ---------------------------------------------------------------------------------------------------------------------------------------"
	if( runMake == 1):
		start_1 = time.time()

		TEMP = 'make -j80 -C '+FARSIGHT_BIN+' > '+LOCAL_DATASET_PATH_LOG +'/runMakeFarsight.log 2>&1'
		TEMP2 = subprocess.Popen(TEMP, shell=True)
		print 'Making Farsight '
		TEMP2.communicate()

		elapsed_1 = (time.time() - start_1)
		print "\t\tTime_1: hm: "+str(round(elapsed_1/3600))+":"+str(round(elapsed_1/60))

	print "# ---------------------------------------------------------------------------------------------------------------------------------------"
	print "# Create Folder: "+DATA_FOLDER+' '+str(datetime.datetime.now())
	print "# ---------------------------------------------------------------------------------------------------------------------------------------"

	# Different logic in case of other structure
	# Different logic in case of other structure
	#if( MOVE_RESULTS == 1 ):
		#if( os.path.exists(GLOBAL_DATASET_PATH_RESULTS) ):
			#shutil.rmtree(GLOBAL_DATASET_PATH_RESULTS)

	a010_CreateFolders.main( LOCAL_DEB_DATASET_PATH, LOCAL_DATASET_PATH_EXE, LOCAL_DATASET_PATH_LOG, LOCAL_DATASET_PATH_DEBUG, LOCAL_DATASET_PATH_DATA, LOCAL_DATASET_PATH_DATA_DEBUG, LOCAL_DATASET_PATH_SEGM, LOCAL_DATASET_PATH_SEGM_DEBUG, LOCAL_DATASET_PATH_SEGM_DEBUG_L2, LOCAL_DATASET_PATH_SEGM_TEMP, LOCAL_DATASET_PATH_TRAC, LOCAL_DATASET_PATH_TRAC_DEBUG, LOCAL_DATASET_PATH_TRAC_DEBUG_L2, LOCAL_DATASET_PATH_TRAC_TEMP, GLOBAL_DATASET_PATH_RESULTS, LOCAL_DATASET_PATH_TRAC_RESULTS, LOCAL_DATASET_PATH_ASTRO_TRAC, LOCAL_DATASET_PATH_ASTRO_TRAC_RESULTS, LOCAL_DATASET_PATH_ASTRO_TRAC_DEBUG, LOCAL_DATASET_PATH_ASTRO_TRAC_DEBUG_L2, LOCAL_DATASET_PATH_ASTRO_TRAC_TEMP)

	if( os.path.exists(GLOBAL_DATASET_PARAMS+'/options_mnt') ):
		shutil.copy(GLOBAL_DATASET_PARAMS+'/options_mnt', LOCAL_DATASET_PATH_PARAMETERS+'/options_mnt')
	if( os.path.exists(GLOBAL_DATASET_PARAMS+'/options_curvelets') ):
		shutil.copy(GLOBAL_DATASET_PARAMS+'/options_curvelets', LOCAL_DATASET_PATH_PARAMETERS+'/options_curvelets')

	print "# ---------------------------------------------------------------------------------------------------------------------------------------"
	print "# Test previous results: "+DATA_FOLDER+' '+str(datetime.datetime.now())
	print "# ---------------------------------------------------------------------------------------------------------------------------------------"
	
	flag = 0
	sameFile = 0
	previousCurve = 0
	if( os.path.exists(GLOBAL_DATASET_PATH_RESULTS+'/Data/options_curvelets') ):
		flag = 1
	if( (flag == 1)):
		if( filecmp.cmp(LOCAL_DATASET_PATH_PARAMETERS+'/options_curvelets', GLOBAL_DATASET_PATH_RESULTS+'/Data/'+'options_curvelets') ):
			sameFile = 1;
	if( sameFile == 1 ):
		print "Curvelets found in: "+GLOBAL_DATASET_PATH_RESULTS+" --> Copy curvelets results"
		TEMP2 = 'scp '+GLOBAL_DATASET_PATH_RESULTS+'/Data/options_curvelets'+' '+LOCAL_DATASET_PATH_DATA+'/options_curvelets'
		print '\t'+TEMP2
		TEMP8 = subprocess.Popen(TEMP2, shell=True)
		TEMP8.communicate()
		TEMP2 = 'scp '+GLOBAL_DATASET_PATH_RESULTS+'/Data/*GFPdsu_BS_CV.nrrd'+' '+LOCAL_DATASET_PATH_DATA+'/'
		print '\t'+TEMP2
		TEMP8 = subprocess.Popen(TEMP2, shell=True)
		TEMP8.communicate()
		FILE_GFP_BS_CV = a019_FindFile.main( LOCAL_DATASET_PATH_DATA, '*GFPdsu_BS_CV.' )
		previousCurve = 1
		runCurv = 0
		

	# Solo si existe
	flag = 0
	sameFile = 0
	if( os.path.exists(GLOBAL_DATASET_PATH_RESULTS+'/Data/options_segmentation') & (previousCurve == 1) ):
		print "Segmentation found in: "+GLOBAL_DATASET_PATH_RESULTS+" --> Copy segmentation results"
		TEMP2 = 'scp '+GLOBAL_DATASET_PATH_RESULTS+'/Data/options_segmentation'+' '+LOCAL_DATASET_PATH_DATA+'/options_segmentation'
		print '\t'+TEMP2
		TEMP8 = subprocess.Popen(TEMP2, shell=True)
		TEMP8.communicate()
		#TEMP2 = 'scp '+GLOBAL_DATASET_PATH_RESULTS+'/Data/*label*'+' '+LOCAL_DATASET_PATH_DATA+'/'
		#print '\t'+TEMP2
		#TEMP8 = subprocess.Popen(TEMP2, shell=True)
		#TEMP8.communicate()
		TEMP2 = 'scp '+GLOBAL_DATASET_PATH_RESULTS+'/Data/*soma*'+' '+LOCAL_DATASET_PATH_DATA+'/'
		print '\t'+TEMP2
		TEMP8 = subprocess.Popen(TEMP2, shell=True)
		TEMP8.communicate()
		#TEMP2 = 'scp '+GLOBAL_DATASET_PATH_RESULTS+'/Data/*GFPdsu_BS_CV_RE.nrrd'+' '+LOCAL_DATASET_PATH_DATA+'/'
		#print '\t'+TEMP2
		#TEMP8 = subprocess.Popen(TEMP2, shell=True)
		#TEMP8.communicate()
		runCopy = 0		# Flag to move images
		runCopy_db = 0
		runMake = 0		# Flag Make Farsight
		runBack = 0		# Flag to run background substraction
		runBack_db = 0

		runCurv = 0		# Flag to run Curvelets
		runCurv_db = 0
		runRescale = 0
		runRescale_db = 0
		runRescale_bit = 0
		runRescale_bit_db = 0
		runSegm = 0		# Flag to run Segmentation
		runSegm_db = 0


	print "# ---------------------------------------------------------------------------------------------------------------------------------------"
	print "# Copy the files to the local dataset: "+DATA_FOLDER+' '+str(datetime.datetime.now())
	print "# ---------------------------------------------------------------------------------------------------------------------------------------"
	if runCopy == 1:
		start_1 = time.time()

		if haveCy5 == 1:
			FILE_Cy5 = a020_CopyFiles.main( GLOBAL_DATASET_PATH, LOCAL_DATASET_PATH_DATA, '*Cy5dsu.' )
		if haveTRT == 1:
			FILE_TRI = a020_CopyFiles.main( GLOBAL_DATASET_PATH, LOCAL_DATASET_PATH_DATA, '*TRITCdsu.' )
		if haveGFP == 1:
			FILE_GFP = a020_CopyFiles.main( GLOBAL_DATASET_PATH, LOCAL_DATASET_PATH_DATA, '*GFPdsu.' )
		if haveDAP == 1:
			FILE_DAP = a020_CopyFiles.main( GLOBAL_DATASET_PATH, LOCAL_DATASET_PATH_DATA, '*DAPIdsu.' )
		elapsed_1 = (time.time() - start_1)
		print "\t\tTime_1: hm: "+str(round(elapsed_1/3600))+":"+str(round(elapsed_1/60))

		if( FILE_Cy5 == '' ):
			haveCy5 = 0
		if( FILE_TRI == '' ):
			haveTRT = 0
		if( FILE_GFP == '' ):
			haveGFP = 0
		if( FILE_DAP == '' ):
			haveDAP = 0

	if runCopy_db == 1:
		start_2 = time.time()
		runCopy_db_log = LOCAL_DATASET_PATH_LOG +'/runCopyProjections_db.log'
		TEMP_FILE = open(runCopy_db_log, 'w')
		TEMP_FILE.write('CopyLog\n')
		TEMP_FILE.close()

		if haveCy5 == 1:
			a021_Project.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_DATA_DEBUG, FILE_Cy5, runCopy_db_log, 'ORG_RES_HISTO', 'TIFF' )
		if haveTRT == 1:
			a021_Project.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_DATA_DEBUG, FILE_TRI, runCopy_db_log, 'ORG_RES_HISTO', 'TIFF' )
		if haveGFP == 1:
			a021_Project.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_DATA_DEBUG, FILE_GFP, runCopy_db_log, 'ORG_RES_HISTO', 'TIFF' )
		if haveDAP == 1:
			a021_Project.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_DATA_DEBUG, FILE_DAP, runCopy_db_log, 'ORG_RES_HISTO', 'TIFF' )
		elapsed_2 = (time.time() - start_2)
		print "\t\tTime_2: hm: "+str(round(elapsed_2/3600))+":"+str(round(elapsed_1/60))

	print "# ---------------------------------------------------------------------------------------------------------------------------------------"
	print "# Run Backgroun Substraction: "+DATA_FOLDER+' '+str(datetime.datetime.now())
	print "# ---------------------------------------------------------------------------------------------------------------------------------------"
	if runBack == 1:
		start_1 = time.time()
		if haveCy5 == 1:
			FILE_Cy5_BS = FILE_Cy5+'_BS'
			a030_RunBack.main( LOCAL_DATASET_PATH_PARAMETERS, FILE_Cy5_BS, FILE_Cy5 )
		if haveTRT == 1:
			FILE_TRI_BS = FILE_TRI+'_BS'
			a030_RunBack.main( LOCAL_DATASET_PATH_PARAMETERS, FILE_TRI_BS, FILE_TRI )
		if haveGFP == 1:
			FILE_GFP_BS = FILE_GFP+'_BS'
			a030_RunBack.main( LOCAL_DATASET_PATH_PARAMETERS, FILE_GFP_BS, FILE_GFP )
		if haveDAP == 1:
			FILE_DAP_BS = FILE_DAP+'_BS'
			a030_RunBack.main( LOCAL_DATASET_PATH_PARAMETERS, FILE_DAP_BS, FILE_DAP )
		elapsed_1 = (time.time() - start_1)
		print "\t\tTime_1: hm: "+str(round(elapsed_1/3600))+":"+str(round(elapsed_1/60))	

	if runBack_db == 1:
		start_2 = time.time()
		runBack_db_log = LOCAL_DATASET_PATH_LOG +'/runBackSubst_db.log'
		TEMP_FILE = open(runBack_db_log, 'w')
		TEMP_FILE.write('BackLog\n')
		TEMP_FILE.close()

		if haveCy5 == 1:
			a021_Project.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_DATA_DEBUG, FILE_Cy5_BS, runBack_db_log, 'ORG_RES_HISTO', 'TIFF' )
		if haveTRT == 1:
			a021_Project.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_DATA_DEBUG, FILE_TRI_BS, runBack_db_log, 'ORG_RES_HISTO', 'TIFF' )
		if haveGFP == 1:
			a021_Project.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_DATA_DEBUG, FILE_GFP_BS, runBack_db_log, 'ORG_RES_HISTO', 'TIFF' )
		if haveDAP == 1:
			a021_Project.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_DATA_DEBUG, FILE_DAP_BS, runBack_db_log, 'ORG_RES_HISTO', 'TIFF' )
		elapsed_2 = (time.time() - start_2)
		print "\t\tTime_2: hm: "+str(round(elapsed_2/3600))+":"+str(round(elapsed_1/60))

	print "# ---------------------------------------------------------------------------------------------------------------------------------------"
	print "# Run Curvelets: "+DATA_FOLDER+' '+str(datetime.datetime.now())
	print "# ---------------------------------------------------------------------------------------------------------------------------------------"
	if runCurv == 1:
		start_1 = time.time()
		FILE_GFP_BS_CV = FILE_GFP_BS+'_CV'
		runCurv_log = LOCAL_DATASET_PATH_LOG +'/runCurvelets.log'
		TEMP_FILE = open(runCurv_log, 'w')
		TEMP_FILE.write('CurvLog\n')
		TEMP_FILE.close()
		a040_RunCurv.main( FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_PARAMETERS, FILE_GFP_BS, runCurv_log, LOCAL_DATASET_PATH_DATA )
		elapsed_1 = (time.time() - start_1)
		print "\t\tTime_1: hm: "+str(round(elapsed_1/3600))+":"+str(round(elapsed_1/60))

	if runCurv_db == 1:
		start_2 = time.time()
		runCurv_db_log = LOCAL_DATASET_PATH_LOG +'/runCurvelets_db.log'
		TEMP_FILE = open(runCurv_db_log, 'w')
		TEMP_FILE.write('CurvLog\n')
		TEMP_FILE.close()
		a021_Project.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_DATA_DEBUG, FILE_GFP_BS_CV, runCurv_db_log, 'ORG_RES_HISTO', 'TIFF' )
		elapsed_2 = (time.time() - start_2)
		print "\t\tTime_2: "+str(elapsed_2/60)

	print "# ---------------------------------------------------------------------------------------------------------------------------------------"
	print "# Run Rescaling: "+DATA_FOLDER+' '+str(datetime.datetime.now())
	print "# ---------------------------------------------------------------------------------------------------------------------------------------"
	if runRescale == 1:
		start_1 = time.time()
		runRescale_log = LOCAL_DATASET_PATH_LOG +'/runRescale.log'
		TEMP_FILE = open(runRescale_log, 'w')
		TEMP_FILE.write('RunRescaleLog\n')
		TEMP_FILE.close()
		if haveCy5 == 1:
			FILE_Cy5_BS_RE = FILE_Cy5_BS+'_RE'
			a050_RunResc.main( FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_PARAMETERS, FILE_Cy5_BS_RE, FILE_Cy5_BS, runRescale_log )
		if haveTRT == 1:
			FILE_TRI_BS_RE = FILE_TRI_BS+'_RE'
			a050_RunResc.main( FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_PARAMETERS, FILE_TRI_BS_RE, FILE_TRI_BS, runRescale_log )
		if haveGFP == 1:
			FILE_GFP_BS_CV_RE = FILE_GFP_BS_CV+'_RE'
			a050_RunResc.main( FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_PARAMETERS, FILE_GFP_BS_CV_RE, FILE_GFP_BS_CV, runRescale_log )
		if haveDAP == 1:
			FILE_DAP_BS_RE = FILE_DAP_BS+'_RE'
			a050_RunResc.main( FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_PARAMETERS, FILE_DAP_BS_RE, FILE_DAP_BS, runRescale_log )
		elapsed_1 = (time.time() - start_1)
		print "\t\tTime_1: hm: "+str(round(elapsed_1/3600))+":"+str(round(elapsed_1/60))	

	if runRescale_db == 1:
		start_2 = time.time()
		runRescale_db_log = LOCAL_DATASET_PATH_LOG +'/runRescale_db.log'
		TEMP_FILE = open(runRescale_db_log, 'w')
		TEMP_FILE.write('RescaleLog\n')
		TEMP_FILE.close()

		if haveCy5 == 1:
			a021_Project.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_DATA_DEBUG, FILE_Cy5_BS_RE, runRescale_db_log, 'ORG_RES_HISTO', 'TIFF' )
		if haveTRT == 1:
			a021_Project.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_DATA_DEBUG, FILE_TRI_BS_RE, runRescale_db_log, 'ORG_RES_HISTO', 'TIFF' )
		if haveGFP == 1:
			a021_Project.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_DATA_DEBUG, FILE_GFP_BS_CV_RE, runRescale_db_log, 'ORG_RES_HISTO', 'TIFF' )
		if haveDAP == 1:
			a021_Project.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_DATA_DEBUG, FILE_DAP_BS_RE, runRescale_db_log, 'ORG_RES_HISTO', 'TIFF' )
		elapsed_2 = (time.time() - start_2)
		print "\t\tTime_2: hm: "+str(round(elapsed_2/3600))+":"+str(round(elapsed_1/60))

	print "# ---------------------------------------------------------------------------------------------------------------------------------------"
	print "# Run Rescaling_bit: "+DATA_FOLDER+' '+str(datetime.datetime.now())
	print "# ---------------------------------------------------------------------------------------------------------------------------------------"
	if runRescale_bit == 1:
		start_1 = time.time()
		runRescale_bit_log = LOCAL_DATASET_PATH_LOG +'/runRescale_bit.log'
		TEMP_FILE = open(runRescale_bit_log, 'w')
		TEMP_FILE.write('RunRescale_bitLog\n')
		TEMP_FILE.close()
		if haveCy5 == 1:
			FILE_Cy5_BS_RE_bit = FILE_Cy5_BS+'_RE_bit'
			a051_RunResc8bit.main( FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_PARAMETERS, FILE_Cy5_BS_RE_bit, FILE_Cy5_BS, runRescale_bit_log )
		if haveTRT == 1:
			FILE_TRI_BS_RE_bit = FILE_TRI_BS+'_RE_bit'
			a051_RunResc8bit.main( FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_PARAMETERS, FILE_TRI_BS_RE_bit, FILE_TRI_BS, runRescale_bit_log )
		if haveGFP == 1:
			FILE_GFP_BS_CV_RE_bit = FILE_GFP_BS_CV+'_RE_bit'
			a051_RunResc8bit.main( FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_PARAMETERS, FILE_GFP_BS_CV_RE_bit, FILE_GFP_BS_CV, runRescale_bit_log )
			a052_RunResc8bitMhd.main( FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_PARAMETERS, FILE_GFP_BS_CV_RE_bit, FILE_GFP_BS_CV, runRescale_bit_log )
		if haveDAP == 1:
			FILE_DAP_BS_RE_bit = FILE_DAP_BS+'_RE_bit'
			a051_RunResc8bit.main( FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_PARAMETERS, FILE_DAP_BS_RE_bit, FILE_DAP_BS, runRescale_bit_log )	
		elapsed_1 = (time.time() - start_1)
		print "\t\tTime_1: hm: "+str(round(elapsed_1/3600))+":"+str(round(elapsed_1/60))

	if runRescale_bit_db == 1:
		start_2 = time.time()
		runRescale_bit_db_log = LOCAL_DATASET_PATH_LOG +'/runRescale_bit_db.log'
		TEMP_FILE = open(runRescale_bit_db_log, 'w')
		TEMP_FILE.write('Rescale_bitLog\n')
		TEMP_FILE.close()

		if haveCy5 == 1:
			a022_Project8bits.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_DATA_DEBUG, FILE_Cy5_BS_RE_bit, runRescale_bit_db_log, 'ORG_RES_HISTO', 'TIFF' )
		if haveTRT == 1:
			a022_Project8bits.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_DATA_DEBUG, FILE_TRI_BS_RE_bit, runRescale_bit_db_log, 'ORG_RES_HISTO', 'TIFF' )
		if haveGFP == 1:
			a022_Project8bits.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_DATA_DEBUG, FILE_GFP_BS_CV_RE_bit, runRescale_bit_db_log, 'ORG_RES_HISTO', 'TIFF' )
		if haveDAP == 1:
			a022_Project8bits.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_DATA_DEBUG, FILE_DAP_BS_RE_bit, runRescale_bit_db_log, 'ORG_RES_HISTO', 'TIFF' )
		elapsed_2 = (time.time() - start_2)
		print "\t\tTime_2: hm: "+str(round(elapsed_2/3600))+":"+str(round(elapsed_1/60))

	print "# ---------------------------------------------------------------------------------------------------------------------------------------"
	print "# Run Segmentation: "+DATA_FOLDER+' '+str(datetime.datetime.now())
	print "# ---------------------------------------------------------------------------------------------------------------------------------------"
	if runSegm == 1:
		start_1 = time.time()
		runSegm_log = LOCAL_DATASET_PATH_LOG +'/runSegmentation.log'
		TEMP_FILE = open(runSegm_log, 'w')
		TEMP_FILE.write('RunSegmentationLog\n')
		TEMP_FILE.close()

		optionsSegm = LOCAL_DATASET_PATH_SEGM +'/options_segmentation'
		TEMP_FILE = open(optionsSegm, 'w')
		if( TEST_RUN == 1):
			TEMP_FILE.write('-xTile 204\n')
			TEMP_FILE.write('-yTile 204\n')
			TEMP_FILE.write('-zTile 102\n')
			TEMP_FILE.write('-xTileBor 20\n')
			TEMP_FILE.write('-yTileBor 20\n')
			TEMP_FILE.write('-zTileBor 10\n')
		else:
			TEMP_FILE.write('-xTile 800\n')
			TEMP_FILE.write('-yTile 800\n')
			TEMP_FILE.write('-zTile 400\n')
			TEMP_FILE.write('-xTileBor 200\n')
			TEMP_FILE.write('-yTileBor 200\n')
			TEMP_FILE.write('-zTileBor 100\n')

		TEMP_FILE.write('-num_threads 80\n')
		if haveCy5 == 1:
			TEMP_FILE.write('-Cy5_Image '+FILE_Cy5_BS_RE_bit+'\n')
		if haveTRT == 1:
			TEMP_FILE.write('-TRI_Image '+FILE_TRI_BS_RE_bit+'\n')
		if haveGFP == 1:
			TEMP_FILE.write('-GFP_Image '+FILE_GFP_BS_CV_RE_bit+'\n')
		if haveDAP == 1:
			TEMP_FILE.write('-DAP_Image '+FILE_DAP_BS_RE_bit+'\n')
		TEMP_FILE.write('-isSmall '+SMALLIMAGE+'\n')
		TEMP_FILE.write('-segParams '+LOCAL_DATASET_PATH_PARAMETERS+'/Seg_Params.ini'+'\n')
		TEMP_FILE.write('-projectDefinition '+LOCAL_DATASET_PATH_PARAMETERS+'/ProjectDefinition.xml'+'\n')
		TEMP_FILE.write('-optionsMNT '+LOCAL_DATASET_PATH_PARAMETERS+'/options_mnt'+'\n')
		TEMP_FILE.write('-outPath '+LOCAL_DATASET_PATH_SEGM+'\n')
		TEMP_FILE.write('-outPathDebug '+LOCAL_DATASET_PATH_SEGM_DEBUG+'\n')
		TEMP_FILE.write('-outPathDebugLevel2 '+LOCAL_DATASET_PATH_SEGM_DEBUG_L2+'\n')
		TEMP_FILE.write('-outPathTemp '+LOCAL_DATASET_PATH_SEGM_TEMP+'\n')
		TEMP_FILE.write('-outPathData '+LOCAL_DATASET_PATH_DATA+'\n')

		TEMP_FILE.close()

		#FILE_Cy5_BS_RE_bit = FILE_Cy5_BS+'_RE_bit'
		a060_RunSegmentation.main( FILE_GFP_BS_CV_RE_bit, FARSIGHT_BIN_EXE, optionsSegm, runSegm_log, LOCAL_DATASET_PATH_DATA )
		elapsed_1 = (time.time() - start_1)
		print "\t\tTime_1: hm: "+str(round(elapsed_1/3600))+":"+str(round(elapsed_1/60))
		
		FILE_LABEL = LOCAL_DATASET_PATH_DATA+'\label'
		FILE_SOMA = LOCAL_DATASET_PATH_DATA+'\soma'
		FILE_SOMA_CEN = LOCAL_DATASET_PATH_DATA+'\soma_centrois.txt'

	if( REMOVE_TEMP_SEGM == 1 ):
		shutil.rmtree(LOCAL_DATASET_PATH_SEGM_TEMP)
		#os.makedirs(LOCAL_DATASET_PATH_SEGM_TEMP)

	if runSegm_db == 1:
		start_2 = time.time()
		runSegm_db_log = LOCAL_DATASET_PATH_LOG +'/runSegm_db.log'
		TEMP_FILE = open(runSegm_db_log, 'w')
		TEMP_FILE.write('SegmentLog\n')
		TEMP_FILE.close()

		a021_Project.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_SEGM_DEBUG, LOCAL_DATASET_PATH_DATA+'/soma', runSegm_db_log, 'ORG_RES_BIN_HISTO', 'TIFF' )
		a021_Project.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_SEGM_DEBUG, LOCAL_DATASET_PATH_DATA+'/label', runSegm_db_log, 'ORG_RES_BIN_HISTO', 'TIFF' )
		a023_ProjectRGB.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_SEGM_DEBUG, FILE_GFP_BS_CV_RE_bit, LOCAL_DATASET_PATH_DATA+'/soma', '_GFP_SOMA_', runSegm_db_log )
		a023_ProjectRGB.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_SEGM_DEBUG, FILE_DAP_BS_RE_bit, LOCAL_DATASET_PATH_DATA+'/label', '_DAPI_LABEL_', runSegm_db_log )
		a023_ProjectRGB.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_SEGM_DEBUG, FILE_GFP_BS_CV_RE_bit, FILE_DAP_BS_RE_bit, 'GF_PDAP_I', runSegm_db_log )

		elapsed_2 = (time.time() - start_2)
		print "\t\tTime_2: hm: "+str(round(elapsed_2/3600))+":"+str(round(elapsed_1/60))
		
	
	print "# ---------------------------------------------------------------------------------------------------------------------------------------"
	print "# Run Tracing: "+DATA_FOLDER+' '+str(datetime.datetime.now())
	print "# ---------------------------------------------------------------------------------------------------------------------------------------"
	if runTrac == 1:
		start_1 = time.time()
		runTrac_log = LOCAL_DATASET_PATH_LOG +'/runTracing.log'
		TEMP_FILE = open(runTrac_log, 'w')
		TEMP_FILE.write('RunTracingLog\n')
		TEMP_FILE.close()

		optionsTracing = LOCAL_DATASET_PATH_TRAC +'/options_tracing'
		TEMP_FILE = open(optionsTracing, 'w')
		if( TEST_RUN == 1):
			TEMP_FILE.write('-xTile 200\n')
			TEMP_FILE.write('-yTile 200\n')
			TEMP_FILE.write('-zTile 50\n')
		else:
			TEMP_FILE.write('-xTile 600\n')
			TEMP_FILE.write('-yTile 600\n')
			TEMP_FILE.write('-zTile 300\n')
		TEMP_FILE.write('-num_threads 80\n')
		if haveCy5 == 1:
			TEMP_FILE.write('-Cy5_Image '+FILE_Cy5_BS_RE_bit+'\n')
		if haveTRT == 1:
			TEMP_FILE.write('-TRI_Image '+FILE_TRI_BS_RE_bit+'\n')
		if haveGFP == 1:
			TEMP_FILE.write('-GFP_Image '+FILE_GFP_BS_CV+'\n')
		if haveDAP == 1:
			TEMP_FILE.write('-DAP_Image '+FILE_DAP_BS_RE_bit+'\n')
		TEMP_FILE.write('-Soma_Centroids '+LOCAL_DATASET_PATH_DATA+'/soma_table_centroids.txt'+'\n')
		TEMP_FILE.write('-Soma_Montage '+LOCAL_DATASET_PATH_DATA+'/soma'+'\n')
		TEMP_FILE.write('-isSmall '+SMALLIMAGE+'\n')
		TEMP_FILE.write('-traceParams '+LOCAL_DATASET_PATH_PARAMETERS+'/options_mnt'+'\n')
		TEMP_FILE.write('-outPath '+LOCAL_DATASET_PATH_TRAC_RESULTS+'\n')
		TEMP_FILE.write('-outPathDebug '+LOCAL_DATASET_PATH_TRAC_DEBUG+'\n')
		TEMP_FILE.write('-outPathDebugLevel2 '+LOCAL_DATASET_PATH_TRAC_DEBUG_L2+'\n')
		TEMP_FILE.write('-outPathTemp '+LOCAL_DATASET_PATH_TRAC_TEMP+'\n')
		TEMP_FILE.close()

		a070_RunTracing.main( FARSIGHT_BIN_EXE, optionsTracing, runTrac_log, LOCAL_DATASET_PATH_DATA )
		elapsed_1 = (time.time() - start_1)
		print "\t\tTime_1: hm: "+str(round(elapsed_1/3600))+":"+str(round(elapsed_1/60))

	if runTrac_db == 1:
		start_2 = time.time()
		runTrac_db_log = LOCAL_DATASET_PATH_LOG +'/runTrac_db.log'
		TEMP_FILE = open(runTrac_db_log, 'w')
		TEMP_FILE.write('TracingLog\n')
		TEMP_FILE.close()

		# Trace editor projection
		a024_ProjectFloat.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_TRAC_DEBUG, LOCAL_DATASET_PATH_TRAC_TEMP+'/GFP_MNT_PRE', runTrac_db_log, 'ORG', 'NRRD' )
		#a021_Project.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_SEGM_DEBUG, FILE_GFP_BS_CV_RE_bit+'_label', runSegm_db_log, 'ORG_RES_BIN' )
		#a023_ProjectRGB.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_SEGM_DEBUG, FILE_GFP_BS_CV_RE_bit, FILE_GFP_BS_CV_RE_bit+'_soma', '_GFP_SOMA_', runSegm_db_log )
		#a023_ProjectRGB.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_SEGM_DEBUG, FILE_DAP_BS_RE_bit, FILE_GFP_BS_CV_RE_bit+'_label', '_DAPI_LABEL_', runSegm_db_log )
		elapsed_2 = (time.time() - start_2)
		print "\t\tTime_2: hm: "+str(round(elapsed_2/3600))+":"+str(round(elapsed_1/60))

	if( REMOVE_TEMP_TRAC == 1 ):
		shutil.rmtree(LOCAL_DATASET_PATH_TRAC_TEMP)

		
	print "# ---------------------------------------------------------------------------------------------------------------------------------------"
	print "# Run AstroTracing: "+DATA_FOLDER+' '+str(datetime.datetime.now())
	print "# ---------------------------------------------------------------------------------------------------------------------------------------"
	if runAstroTrac == 1:
		start_1 = time.time()
		runAstroTrac_log = LOCAL_DATASET_PATH_LOG +'/runAstroTracing.log'
		TEMP_FILE = open(runAstroTrac_log, 'w')
		TEMP_FILE.write('runAstroTracing\n')
		TEMP_FILE.close()
		
		runDistMap_log = LOCAL_DATASET_PATH_LOG +'/runDistMap.log'
		TEMP_FILE = open(runDistMap_log, 'w')
		TEMP_FILE.write('runDistMap\n')
		TEMP_FILE.close()
		
		FILE_LABEL_DIST_MAP = LOCAL_DATASET_PATH_DATA+'/label_dist_map'
		a061_runDistMap.main( FARSIGHT_BIN_EXE, FILE_LABEL_DIST_MAP, FILE_LABEL, runDistMap_log )
		
		optionsAstroTracing = LOCAL_DATASET_PATH_ASTRO_TRAC +'/options_astrocyte_tracing'
			TEMP_FILE = open(optionsAstroTracing, 'w')
			if( TEST_RUN == 1):
				TEMP_FILE.write('-xTile 204\n')
				TEMP_FILE.write('-yTile 204\n')
				TEMP_FILE.write('-zTile 102\n')
				TEMP_FILE.write('-xTileBor 20\n')
				TEMP_FILE.write('-yTileBor 20\n')
				TEMP_FILE.write('-zTileBor 10\n')
			else:
				TEMP_FILE.write('-xTile 800\n')
				TEMP_FILE.write('-yTile 800\n')
				TEMP_FILE.write('-zTile 400\n')
				TEMP_FILE.write('-xTileBor 200\n')
				TEMP_FILE.write('-yTileBor 200\n')
				TEMP_FILE.write('-zTileBor 100\n')
			TEMP_FILE.write('-num_threads 80\n')
			if haveCy5 == 1:
				TEMP_FILE.write('-Cy5_Image '+FILE_Cy5_BS_RE_bit+'\n')
			if haveTRT == 1:
				TEMP_FILE.write('-TRI_Image '+FILE_TRI_BS_RE_bit+'\n')
			if haveGFP == 1:
				TEMP_FILE.write('-GFP_Image '+FILE_GFP_BS_CV+'\n')
			if haveDAP == 1:
				TEMP_FILE.write('-DAP_Image '+FILE_DAP_BS_RE_bit+'\n')
			TEMP_FILE.write('-Dist_Map_Image '+FILE_LABEL_DIST_MAP+'\n')
			TEMP_FILE.write('-Soma_Centroids '+LOCAL_DATASET_PATH_DATA+'/soma_table_centroids.txt'+'\n')
			TEMP_FILE.write('-Soma_Montage '+LOCAL_DATASET_PATH_DATA+'/soma'+'\n')
			TEMP_FILE.write('-isSmall '+SMALLIMAGE+'\n')
			TEMP_FILE.write('-astroTraceParams '+LOCAL_DATASET_PATH_PARAMETERS+'/options_as_tr'+'\n')
			TEMP_FILE.write('-outPath '+LOCAL_DATASET_PATH_ASTRO_TRAC_RESULTS+'\n')
			TEMP_FILE.write('-outPathDebug '+LOCAL_DATASET_PATH_ASTRO_TRAC_DEBUG+'\n')
			TEMP_FILE.write('-outPathDebugLevel2 '+LOCAL_DATASET_PATH_ASTRO_TRAC_DEBUG_L2+'\n')
			TEMP_FILE.write('-outPathTemp '+LOCAL_DATASET_PATH_ASTRO_TRAC_TEMP+'\n')
		TEMP_FILE.close()

		a075_RunAstroTracing.main( FARSIGHT_BIN_EXE, optionsAstroTracing, runAstroTrac_log, LOCAL_DATASET_PATH_DATA )
		elapsed_1 = (time.time() - start_1)
		print "\t\tTime_1: hm: "+str(round(elapsed_1/3600))+":"+str(round(elapsed_1/60))

	if runAstroTrac_db == 1:
		start_2 = time.time()
		runAstroTrac_db_log = LOCAL_DATASET_PATH_LOG +'/runAstroTrac_db.log'
		TEMP_FILE = open(runAstroTrac_db_log, 'w')
		TEMP_FILE.write('AstroTracingLog\n')
		TEMP_FILE.close()

		# Trace editor projection
		#a024_ProjectFloat.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_ASTRO_TRAC_DEBUG, LOCAL_DATASET_PATH_ASTRO_TRAC_TEMP+'/GFP_MNT_PRE', runAstroTrac_db_log, 'ORG', 'NRRD' )
		#a021_Project.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_SEGM_DEBUG, FILE_GFP_BS_CV_RE_bit+'_label', runSegm_db_log, 'ORG_RES_BIN' )
		#a023_ProjectRGB.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_SEGM_DEBUG, FILE_GFP_BS_CV_RE_bit, FILE_GFP_BS_CV_RE_bit+'_soma', '_GFP_SOMA_', runSegm_db_log )
		#a023_ProjectRGB.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_SEGM_DEBUG, FILE_DAP_BS_RE_bit, FILE_GFP_BS_CV_RE_bit+'_label', '_DAPI_LABEL_', runSegm_db_log )
		elapsed_2 = (time.time() - start_2)
		print "\t\tTime_2: hm: "+str(round(elapsed_2/3600))+":"+str(round(elapsed_1/60))
		
		
		
		
		
		
		
	print "# ---------------------------------------------------------------------------------------------------------------------------------------"
	print "# Move Files: "+DATA_FOLDER+' '+str(datetime.datetime.now())
	print "# ---------------------------------------------------------------------------------------------------------------------------------------"


		#if( REMOVE_TEMP_SEGM == 1 ):
			#shutil.rmtree(LOCAL_DATASET_PATH_SEGM_TEMP)
		#if( REMOVE_TEMP_TRAC == 1 ):
			#shutil.rmtree(LOCAL_DATASET_PATH_TRAC_TEMP)

			#shutil.rmtree(LOCAL_DATASET_PATH_SEGM_TEMP)

		#shutil.copytree(LOCAL_DATASET_PATH,GLOBAL_DATASET_PATH_RESULTS)

	# ----------------------------------------------

		#print "Copying files"
		#os.makedirs(GLOBAL_DATASET_PATH_RESULTS+'/Data')
		#dirList=os.listdir(LOCAL_DATASET_PATH_DATA)
		#for file in dirList:
			##f = open(LOCAL_DATASET_PATH_DATA+'/'+file,'r')
			##f.close()
			#shutil.copy(LOCAL_DATASET_PATH_DATA+'/'+file, GLOBAL_DATASET_PATH_RESULTS+'/Data/'+file)

		#os.makedirs(GLOBAL_DATASET_PATH_RESULTS+'/Exe')
		#dirList=os.listdir(LOCAL_DATASET_PATH_EXE)
		#for file in dirList:
			##f = open(LOCAL_DATASET_PATH_EXE+'/'+file,'r')
			##f.close()
			#shutil.copy(LOCAL_DATASET_PATH_EXE+'/'+file, GLOBAL_DATASET_PATH_RESULTS+'/Exe/'+file)

		#os.makedirs(GLOBAL_DATASET_PATH_RESULTS+'/Log')
		#dirList=os.listdir(LOCAL_DATASET_PATH_LOG)
		#for file in dirList:
			##f = open(LOCAL_DATASET_PATH_LOG+'/'+file,'r')
			##f.close()
			#shutil.copy(LOCAL_DATASET_PATH_LOG+'/'+file, GLOBAL_DATASET_PATH_RESULTS+'/Log/'+file)

		#os.makedirs(GLOBAL_DATASET_PATH_RESULTS+'/Parameters')
		#dirList=os.listdir(LOCAL_DATASET_PATH_PARAMETERS)
		#for file in dirList:
			##f = open(LOCAL_DATASET_PATH_PARAMETERS+'/'+file,'r')
			##f.close()
			#shutil.copy(LOCAL_DATASET_PATH_PARAMETERS+'/'+file, GLOBAL_DATASET_PATH_RESULTS+'/Parameters/'+file)

		#os.makedirs(GLOBAL_DATASET_PATH_RESULTS+'/Segm')
		#dirList=os.listdir(LOCAL_DATASET_PATH_SEGM)
		#for file in dirList:
			##f = open(LOCAL_DATASET_PATH_SEGM+'/'+file,'r')
			##f.close()
			#shutil.copy(LOCAL_DATASET_PATH_SEGM+'/'+file, GLOBAL_DATASET_PATH_RESULTS+'/Segm/'+file)

		#os.makedirs(GLOBAL_DATASET_PATH_RESULTS+'/Trace/Results')
		#dirList=os.listdir(LOCAL_DATASET_PATH_TRAC_RESULTS)
		#for file in dirList:
			##f = open(LOCAL_DATASET_PATH_TRAC+'/'+file,'r')
			##f.close()
			#shutil.copy(LOCAL_DATASET_PATH_TRAC_RESULTS+'/'+file, GLOBAL_DATASET_PATH_RESULTS+'/Trace/Results/'+file)

		#shutil.copy(LOCAL_DATASET_PATH_TRAC +'/options_tracing', GLOBAL_DATASET_PATH_RESULTS+'/Trace/options_tracing')

	if( (MOVE_RESULTS == 1) or (MOVE_LOCALLY == 1) ):
		if( MOVE_INTER_STEPS == 0):
			print "Removing data files"

			if haveCy5 == 1:
				if( os.path.exists(FILE_Cy5+'.nrrd') ):
					os.remove(FILE_Cy5+'.nrrd')
				if( os.path.exists(FILE_Cy5_BS+'.nrrd') ):
					os.remove(FILE_Cy5_BS+'.nrrd')
				if( os.path.exists(FILE_Cy5_BS_RE+'.nrrd') ):
					os.remove(FILE_Cy5_BS_RE+'.nrrd')
				if( os.path.exists(FILE_Cy5_BS_RE_bit+'.nrrd') ):
					os.remove(FILE_Cy5_BS_RE_bit+'.nrrd')
			if haveTRT == 1:
				if( os.path.exists(FILE_TRI+'.nrrd') ):
					os.remove(FILE_TRI+'.nrrd')
				if( os.path.exists(FILE_TRI_BS+'.nrrd') ):
					os.remove(FILE_TRI_BS+'.nrrd')
				if( os.path.exists(FILE_TRI_BS_RE+'.nrrd') ):
					os.remove(FILE_TRI_BS_RE+'.nrrd')
				if( os.path.exists(FILE_TRI_BS_RE_bit+'.nrrd') ):
					os.remove(FILE_TRI_BS_RE_bit+'.nrrd')
			if haveGFP == 1:
				if( os.path.exists(FILE_GFP+'.nrrd') ):
					os.remove(FILE_GFP+'.nrrd')
				if( os.path.exists(FILE_GFP_BS+'.nrrd') ):
					os.remove(FILE_GFP_BS+'.nrrd')
				#if( os.path.exists(FILE_GFP_BS_CV+'.nrrd') ):
					#os.remove(FILE_GFP_BS_CV+'.nrrd')
				if( os.path.exists(FILE_GFP_BS_CV_RE+'.nrrd') ):
					os.remove(FILE_GFP_BS_CV_RE+'.nrrd')
				if( os.path.exists(FILE_GFP_BS_CV_RE_bit+'.nrrd') ):
					os.remove(FILE_GFP_BS_CV_RE_bit+'.nrrd')
			if haveDAP == 1:
				if( os.path.exists(FILE_DAP+'.nrrd') ):
					os.remove(FILE_DAP+'.nrrd')
				if( os.path.exists(FILE_DAP_BS+'.nrrd') ):
					os.remove(FILE_DAP_BS+'.nrrd')
				if( os.path.exists(FILE_DAP_BS_RE+'.nrrd') ):
					os.remove(FILE_DAP_BS_RE+'.nrrd')
				if( os.path.exists(FILE_DAP_BS_RE_bit+'.nrrd') ):
					os.remove(FILE_DAP_BS_RE_bit+'.nrrd')

	# ----------------------------------------------
	if( MOVE_RESULTS == 1):
		print "# ---------------------------------------------------------------------------------------------------------------------------------------"
		print "# Check Sums Files and Send: "+DATA_FOLDER+' '+str(datetime.datetime.now())
		print "# ---------------------------------------------------------------------------------------------------------------------------------------"
		a080_RunCheckSumNew.main( LOCAL_DATASET_PATH_DATA, REMOVE_MONTAGES, GLOBAL_DATASET_PATH_RESULTS )

		print "# ---------------------------------------------------------------------------------------------------------------------------------------"
		print "# Move Files to FSdata: "+DATA_FOLDER+' '+str(datetime.datetime.now())
		print "# ---------------------------------------------------------------------------------------------------------------------------------------"
		
		#STARTOFCOPY = str(datetime.datetime.now())
		#print "Copy to global results"
		#if( REMOVE_MONTAGES == 1 ):
			#TEMP2 = 'scp -r '+LOCAL_DATASET_PATH_DATA+' '+GLOBAL_DATASET_PATH_RESULTS+'/;rm -r '+LOCAL_DATASET_PATH_DATA
			#print '\t'
			#TEMP2 = TEMP2 #+';echo "------------------->>>>>>>>>>>>>>>>>>> FINCOPYINGDATA <<<<<<<<<<<<<<<<<<<<-----------------"OF: '+DATA_FOLDER+", START: "+STARTOFCOPY
			#TEMP3 = subprocess.Popen(TEMP2, shell=True)
		#else:
			#TEMP2 = 'scp -r '+LOCAL_DATASET_PATH_DATA+' '+GLOBAL_DATASET_PATH_RESULTS+'/'
			#print '\t'+TEMP2
			#TEMP3 = subprocess.Popen(TEMP2, shell=True)
		#subprocess.call(['scp -r '])
		#TEMP3.communicate()
		TEMP2 = 'scp -r '+LOCAL_DATASET_PATH_EXE+' far-01:'+GLOBAL_DATASET_PATH_RESULTS+'/'
		print '\t'+TEMP2
		TEMP4 = subprocess.Popen(TEMP2, shell=True)
		#TEMP4.communicate()
		TEMP2 = 'scp -r '+LOCAL_DATASET_PATH_LOG+' far-01:'+GLOBAL_DATASET_PATH_RESULTS+'/'
		print '\t'+TEMP2
		TEMP5 = subprocess.Popen(TEMP2, shell=True)
		#TEMP5.communicate()
		TEMP2 = 'scp -r '+LOCAL_DATASET_PATH_PARAMETERS+' far-01:'+GLOBAL_DATASET_PATH_RESULTS+'/'
		print '\t'+TEMP2
		TEMP6 = subprocess.Popen(TEMP2, shell=True)
		#TEMP6.communicate()
		TEMP2 = 'scp -r '+LOCAL_DATASET_PATH_SEGM+' far-01:'+GLOBAL_DATASET_PATH_RESULTS+'/'
		print '\t'+TEMP2
		TEMP7 = subprocess.Popen(TEMP2, shell=True)
		#TEMP7.communicate()
		TEMP2 = 'scp -r '+LOCAL_DATASET_PATH_TRAC+' far-01:'+GLOBAL_DATASET_PATH_RESULTS+'/'
		print '\t'+TEMP2
		TEMP8 = subprocess.Popen(TEMP2, shell=True)
		#TEMP8.communicate()


		#shutil.copytree(LOCAL_DATASET_PATH_DATA,GLOBAL_DATASET_PATH_RESULTS+'/Data') #// Genera errores
		#shutil.copytree(LOCAL_DATASET_PATH_EXE,GLOBAL_DATASET_PATH_RESULTS+'/Exe')
		#shutil.copytree(LOCAL_DATASET_PATH_LOG,GLOBAL_DATASET_PATH_RESULTS+'/Log')
		#shutil.copytree(LOCAL_DATASET_PATH_PARAMETERS,GLOBAL_DATASET_PATH_RESULTS+'/Parameters')
		#shutil.copytree(LOCAL_DATASET_PATH_SEGM,GLOBAL_DATASET_PATH_RESULTS+'/Segm')
		#shutil.copytree(LOCAL_DATASET_PATH_TRACE_SOMAS,GLOBAL_DATASET_PATH_RESULTS+'/TracesAndSomas')

	#if( REMOVE_MONTAGES == 1):
		##shutil.rmtree(LOCAL_DATASET_PATH_DATA)
		
		#if( MOVE_INTER_STEPS == 0):
			#print "Removign inter steps, and "
			#if haveCy5 == 1:
				#os.remove(FILE_Cy5+'.nrrd')
				#os.remove(FILE_Cy5_BS+'.nrrd')
				#os.remove(FILE_Cy5_BS_RE+'.nrrd')
				#os.remove(FILE_Cy5_BS_RE_bit+'.nrrd')
			#if haveTRT == 1:
				#os.remove(FILE_TRI+'.nrrd')
				#os.remove(FILE_TRI_BS+'.nrrd')
				#os.remove(FILE_TRI_BS_RE+'.nrrd')
				#os.remove(FILE_TRI_BS_RE_bit+'.nrrd')
			#if haveGFP == 1:
				#os.remove(FILE_GFP+'.nrrd')
				#os.remove(FILE_GFP_BS+'.nrrd')
				##os.remove(FILE_GFP_BS_CV+'.nrrd')
				#os.remove(FILE_GFP_BS_CV_RE+'.nrrd')
				##os.remove(FILE_GFP_BS_CV_RE_bit+'.nrrd')
			#if haveDAP == 1:
				#os.remove(FILE_DAP+'.nrrd')
				#os.remove(FILE_DAP_BS+'.nrrd')
				#os.remove(FILE_DAP_BS_RE+'.nrrd')
				#os.remove(FILE_DAP_BS_RE_bit+'.nrrd')
	
	#if( REMOVE_NRRD == 1):
		#dirListGlo=os.listdir(LOCAL_DEB_DATASET_PATH)
		#for fileGlo in dirListGlo:
			#FILE = LOCAL_DATASET_PATH_DATA+'/'+fileGlo
			#print "Erasing: "+FILE
			#if fnmatch.fnmatch(fileGlo, '*.nrrd'):
				#os.remove(FILE)

	print "# ---------------------------------------------------------------------------------------------------------------------------------------"
	print "# Move internally: "+DATA_FOLDER+' '+str(datetime.datetime.now())
	print "# ---------------------------------------------------------------------------------------------------------------------------------------"
	if( MOVE_LOCALLY == 1):
		#TEMP2 = 'rm -r '+LOCAL_DATASET_PATH_DATA
		#print '\t'+TEMP2
		#TEMP9 = subprocess.Popen(TEMP2, shell=True)
		#TEMP9.communicate()
		#TEMP2 = 'scp -r '+LOCAL_DATASET_PATH+'/* '+LOCAL_DEB_DATASET_PATH+'/'
		#print '\t'+TEMP2
		#TEMP9 = subprocess.Popen(TEMP2, shell=True)
		#TEMP9.communicate()
		#TEMP2 = 'rm -r '+LOCAL_DEB_DATASET_PATH+'/Data'
		#print '\t'+TEMP2
		#TEMP9 = subprocess.Popen(TEMP2, shell=True)
		#TEMP9.communicate()
		#TEMP2 = 'rm -r '+LOCAL_DATASET_PATH_DEBUG+'/Data'
		#print '\t'+TEMP2
		#TEMP9 = subprocess.Popen(TEMP2, shell=True)
		#TEMP9.communicate()

		TEMP2 = 'scp -r '+LOCAL_DATASET_PATH_DEBUG+' '+LOCAL_DEB_DATASET_PATH+'/'
		print '\t'+TEMP2
		TEMP9 = subprocess.Popen(TEMP2, shell=True)
		TEMP9.communicate()
		TEMP2 = 'scp -r '+LOCAL_DATASET_PATH_EXE+' '+LOCAL_DEB_DATASET_PATH+'/'
		print '\t'+TEMP2
		TEMP10 = subprocess.Popen(TEMP2, shell=True)
		TEMP10.communicate()
		TEMP2 = 'scp -r '+LOCAL_DATASET_PATH_LOG+' '+LOCAL_DEB_DATASET_PATH+'/'
		print '\t'+TEMP2
		TEMP11 = subprocess.Popen(TEMP2, shell=True)
		TEMP11.communicate()
		TEMP2 = 'scp -r '+LOCAL_DATASET_PATH_PARAMETERS+' '+LOCAL_DEB_DATASET_PATH+'/'
		print '\t'+TEMP2
		TEMP12 = subprocess.Popen(TEMP2, shell=True)
		TEMP12.communicate()
		TEMP2 = 'scp -r '+LOCAL_DATASET_PATH_SEGM+' '+LOCAL_DEB_DATASET_PATH+'/'
		print '\t'+TEMP2
		TEMP13 = subprocess.Popen(TEMP2, shell=True)
		TEMP13.communicate()
		TEMP2 = 'scp -r '+LOCAL_DATASET_PATH_TRAC+' '+LOCAL_DEB_DATASET_PATH+'/'
		print '\t'+TEMP2
		TEMP14 = subprocess.Popen(TEMP2, shell=True)
		TEMP14.communicate()

		print "---- >>> DONE DONE <<< ---"

		# NOT REMOVE THE DEBUG FOR NOW
		#TEMP2 = 'rm -r '+LOCAL_DATASET_PATH_DEBUG
		#print '\t'+TEMP2
		#TEMP15 = subprocess.Popen(TEMP2, shell=True)
		#TEMP15.communicate()

