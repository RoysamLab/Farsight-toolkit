#!/usr/bin/env python
# -*- coding: utf-8 -*-

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


import a010_CreateFolders
import a020_CopyFiles
import a021_Project
import a022_Project8bits
import a023_ProjectRGB
import a030_RunBack
import a040_RunCurv
import a050_RunResc
import a051_RunResc8bit
import a060_RunSegmentation


class MyWriter:
	def __init__(self,stdout, filename):
		self.stdout = stdout
		self.logfile = file(filename, 'w')
	def write(self, text):
		self.stdout.write(text)
		self.logfile.write(text)

	def close(self):
		self.stdout.close()
		self.logfile.close()

writerLog = MyWriter(sys.stdout, 'logFile_3.log')
sys.stdout = writerLog




#def main(DATA_FOLDER):

# Parameters
TRY = 01
#DEBUG = 1
#----------------------------------
REMOVE_MONTAGES = 1	# This flag is set in case we want the montages to be removed after the process is done, especially when running many montages in serial we want to make sure not to
MOVE_RESULTS = 0	# If 0 the results will be keep
			# if 1 the results will be moved
			# if 2 the results will be copied (keep and copy to FSDATA)
			# if 3 move everysingle file, exept the folder, which are copied
MOVE_LOCALLY = 1
REMOVE_TEMP_SEG = 1
#----------------------------------
MOVE_INTER_STEPS = 0

SMALLIMAGE = '1'	# if the image is small

runCopy = 1		# Flag to move images
runCopy_db = 1
runMake = 1		# Flag Make Farsight
runBack = 1		# Flag to run background substraction
runBack_db = 1
runCurv = 0#1		# Flag to run Curvelets
runCurv_db = 0#1
runRescale = 0#1
runRescale_db = 0#1
runRescale_bit = 0#1
runRescale_bit_db = 0#1
runSegm = 0#1		# Flag to run Segmentation
runSegm_db = 0#1
runTrac = 0#1		# Flag to run Tracing
runTrac_db = 0#1

haveCy5 = 0
haveTRT = 0
haveGFP = 1
haveDAP = 1

# ---------------------------------------------------------------------------------------------------------------------------------------

SERVER = 'none'
p = os.uname()
if( p[1] == 'Farsight-05.EE.UH.EDU' ):
	SERVER = 'far05'
if( p[1] == 'Farsight-04.EE.UH.EDU' ):
	SERVER = 'far04'

#DATA_FOLDER_ALL = ['/0131_test','/0131_test2'] # For testing dont forget the xTile params
#DATA_FOLDER_ALL = ['/0131_test'] # For testing dont forget the xTile params
#DATA_FOLDER_ALL = ['/0113_NRRD_CROPPED']

#DATA_FOLDER_ALL = ['/0113_NRRD_CROPPED','/0117_NRRD_CROPPED','/0120_NRRD_CROPPED','/0123_NRRD_CROPPED','/0128_NRRD_CROPPED','/0131_NRRD_CROPPED','/0323_NRRD_CROPPED','/0405_NRRD_CROPPED','/0409_NRRD_CROPPED','/0410_NRRD_CROPPED','/0412_NRRD_CROPPED','/1206_NRRD_CROPPED']
#DATA_FOLDER_ALL = ['/0113_NRRD_CROPPED','/0117_NRRD_CROPPED','/0120_NRRD_CROPPED','/0123_NRRD_CROPPED','/0128_NRRD_CROPPED','/0131_NRRD_CROPPED']
#DATA_FOLDER_ALL = ['/0323_NRRD_CROPPED','/0405_NRRD_CROPPED','/0409_NRRD_CROPPED','/0410_NRRD_CROPPED','/0412_NRRD_CROPPED','/1206_NRRD_CROPPED']

#DATA_FOLDER_ALL = ['/0113_NRRD','/0117_NRRD','/0120_NRRD','/0123_NRRD','/0128_NRRD','/0131_NRRD','/0323_NRRD','/0405_NRRD','/0409_NRRD','/0410_NRRD','/0412_NRRD','/1206_NRRD']

DATA_FOLDER_ALL = ['/0131_test']
if( SERVER == 'far04' ):
	DATA_FOLDER_ALL = ['/0113_NRRD','/0117_NRRD','/0120_NRRD','/0123_NRRD','/0128_NRRD','/0131_NRRD','/0323_NRRD','/0405_NRRD','/0409_NRRD','/0410_NRRD','/0412_NRRD','/1206_NRRD']
	#DATA_FOLDER_ALL = ['/0113_NRRD','/0117_NRRD','/0120_NRRD','/0123_NRRD','/0128_NRRD','/0131_NRRD']
	#DATA_FOLDER_ALL = ['/0113_NRRD']
	#DATA_FOLDER_ALL = ['/0131_test']
	#DATA_FOLDER_ALL = ['/0131_test2']
	#DATA_FOLDER_ALL = ['/0131_test3']
	#DATA_FOLDER_ALL = ['/0131_test4']
if( SERVER == 'far05' ):
	#DATA_FOLDER_ALL = ['/0323_NRRD','/0405_NRRD','/0409_NRRD','/0410_NRRD','/0412_NRRD','/1206_NRRD']
	#DATA_FOLDER_ALL = ['/0131_test']
	#DATA_FOLDER_ALL = ['/0131_test2']
	#DATA_FOLDER_ALL = ['/0131_test3']
	DATA_FOLDER_ALL = ['/0131_test4']

# ---------------------------------------------------------------------------------------------------------------------------------------

#def main(DATA_FOLDER):

for DATA_FOLDER in DATA_FOLDER_ALL:

	# ---------------------------------------------------------------------------------------------------------------------------------------
	# Folders names 
	# ---------------------------------------------------------------------------------------------------------------------------------------

	FARSIGHT_BIN = '/data/nicolas/farsight_updated/bin'
	FARSIGHT_BIN_EXE = '/data/nicolas/farsight_updated/bin/exe'

	MAIN_DATA_FOLDER = '/data/nicolas/dataNew'
	MAIN_DEB_DATA_FOLDER = '/data/nicolas/deb'

	LOCAL_DATASET_PATH = MAIN_DATA_FOLDER+DATA_FOLDER
	LOCAL_DEB_DATASET_PATH = MAIN_DEB_DATA_FOLDER+DATA_FOLDER+'_RESULTS_V3_'+SERVER+'_'+str(TRY)

	LOCAL_DATASET_PATH_PARAMETERS = LOCAL_DATASET_PATH+'/Parameters'	# ---> This directory has to exist
	LOCAL_DATASET_PATH_EXE = LOCAL_DATASET_PATH+'/Exe'
	LOCAL_DATASET_PATH_LOG = LOCAL_DATASET_PATH+'/Log'
	LOCAL_DATASET_PATH_DEBUG = LOCAL_DATASET_PATH+'/DEBUG'
	LOCAL_DATASET_PATH_TRACE_SOMAS = LOCAL_DATASET_PATH+'/TracesAndSomas'
	#LOCAL_DATASET_PATH_TRACE_SOMASDDIVIDED = LOCAL_DATASET_PATH+'/TracesAndSomasDivided'
	LOCAL_DATASET_PATH_DATA = LOCAL_DATASET_PATH+'/Data'
	LOCAL_DATASET_PATH_DATA_DEBUG = LOCAL_DATASET_PATH_DEBUG+'/Data'
	LOCAL_DATASET_PATH_SEGM = LOCAL_DATASET_PATH+'/Segm'
	LOCAL_DATASET_PATH_SEGM_DEBUG = LOCAL_DATASET_PATH_DEBUG+'/Segm'
	LOCAL_DATASET_PATH_SEGM_DEBUG_L2 = LOCAL_DATASET_PATH_SEGM_DEBUG+'/Level2'
	LOCAL_DATASET_PATH_SEGM_TEMP = LOCAL_DATASET_PATH_SEGM+'/Temp'

	GLOBAL_DATASET_PATH = "/FSdata.auto/data/DARPA_MOSAICS"+DATA_FOLDER
	#GLOBAL_DATASET_PATH = "/FSdata.auto/data"+DATA_FOLDER
	GLOBAL_DATASET_PATH_RESULTS = '/FSdata.auto/data/DARPA_RESULTS'+DATA_FOLDER+'_RESULTS_V3_'+SERVER+'_'+str(TRY)

	print "# ---------------------------------------------------------------------------------------------------------------------------------------"
	print " Make Farsight: "+DATA_FOLDER
	print "# ---------------------------------------------------------------------------------------------------------------------------------------"
	start_1 = time.clock()

	TEMP = 'make -j80 -C '+FARSIGHT_BIN+' > '+LOCAL_DATASET_PATH_LOG +'/runMakeFarsight.log 2>&1'
	TEMP2 = subprocess.Popen(TEMP, shell=True)
	print 'Making Farsight '
	TEMP2.communicate()

	elapsed_1 = (time.clock() - start_1)
	print "\t\tTime_1: "+str(elapsed_1)

	print "# ---------------------------------------------------------------------------------------------------------------------------------------"
	print "# Create Folder: "+DATA_FOLDER
	print "# ---------------------------------------------------------------------------------------------------------------------------------------"

	# Different logic in case of other structure
	# Different logic in case of other structure
	if( MOVE_RESULTS == 1 ):
		if( os.path.exists(GLOBAL_DATASET_PATH_RESULTS) ):
			shutil.rmtree(GLOBAL_DATASET_PATH_RESULTS)

	a010_CreateFolders.main( LOCAL_DEB_DATASET_PATH, LOCAL_DATASET_PATH_EXE, LOCAL_DATASET_PATH_LOG, LOCAL_DATASET_PATH_DEBUG, LOCAL_DATASET_PATH_TRACE_SOMAS, LOCAL_DATASET_PATH_DATA, LOCAL_DATASET_PATH_DATA_DEBUG, LOCAL_DATASET_PATH_SEGM, LOCAL_DATASET_PATH_SEGM_DEBUG, LOCAL_DATASET_PATH_SEGM_DEBUG_L2, LOCAL_DATASET_PATH_SEGM_TEMP, GLOBAL_DATASET_PATH_RESULTS )

	print "# ---------------------------------------------------------------------------------------------------------------------------------------"
	print "# Copy the files to the local dataset: "+DATA_FOLDER
	print "# ---------------------------------------------------------------------------------------------------------------------------------------"
	if runCopy == 1:
		start_1 = time.clock()

		if haveCy5 == 1:
			FILE_Cy5 = a020_CopyFiles.main( GLOBAL_DATASET_PATH, LOCAL_DATASET_PATH_DATA, '*Cy5dsu.' )
		if haveTRT == 1:
			FILE_TRI = a020_CopyFiles.main( GLOBAL_DATASET_PATH, LOCAL_DATASET_PATH_DATA, '*TRITCdsu.' )
		if haveGFP == 1:
			FILE_GFP = a020_CopyFiles.main( GLOBAL_DATASET_PATH, LOCAL_DATASET_PATH_DATA, '*GFPdsu.' )
		if haveDAP == 1:
			FILE_DAP = a020_CopyFiles.main( GLOBAL_DATASET_PATH, LOCAL_DATASET_PATH_DATA, '*DAPIdsu.' )
		elapsed_1 = (time.clock() - start_1)
		print "\t\tTime_1: "+str(elapsed_1)

	if runCopy_db == 1:
		start_2 = time.clock()
		runCopy_db_log = LOCAL_DATASET_PATH_LOG +'/runCopyProjections_db.log'
		TEMP_FILE = open(runCopy_db_log, 'w')
		TEMP_FILE.write('CopyLog\n')
		TEMP_FILE.close()

		if haveCy5 == 1:
			a021_Project.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_DATA_DEBUG, FILE_Cy5, runCopy_db_log, 'ORG_RES' )
		if haveTRT == 1:
			a021_Project.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_DATA_DEBUG, FILE_TRI, runCopy_db_log, 'ORG_RES' )
		if haveGFP == 1:
			a021_Project.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_DATA_DEBUG, FILE_GFP, runCopy_db_log, 'ORG_RES' )
		if haveDAP == 1:
			a021_Project.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_DATA_DEBUG, FILE_DAP, runCopy_db_log, 'ORG_RES' )
		elapsed_2 = (time.clock() - start_2)
		print "\t\tTime_2: "+str(elapsed_2)

	print "# ---------------------------------------------------------------------------------------------------------------------------------------"
	print "# Run Backgroun Substraction: "+DATA_FOLDER
	print "# ---------------------------------------------------------------------------------------------------------------------------------------"
	if runBack == 1:
		start_1 = time.clock()
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
		elapsed_1 = (time.clock() - start_1)
		print "\t\tTime_1: "+str(elapsed_1)	

	if runBack_db == 1:
		start_2 = time.clock()
		runBack_db_log = LOCAL_DATASET_PATH_LOG +'/runBackSubst_db.log'
		TEMP_FILE = open(runBack_db_log, 'w')
		TEMP_FILE.write('BackLog\n')
		TEMP_FILE.close()

		if haveCy5 == 1:
			a021_Project.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_DATA_DEBUG, FILE_Cy5_BS, runBack_db_log, 'ORG_RES' )
		if haveTRT == 1:
			a021_Project.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_DATA_DEBUG, FILE_TRI_BS, runBack_db_log, 'ORG_RES' )
		if haveGFP == 1:
			a021_Project.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_DATA_DEBUG, FILE_GFP_BS, runBack_db_log, 'ORG_RES' )
		if haveDAP == 1:
			a021_Project.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_DATA_DEBUG, FILE_DAP_BS, runBack_db_log, 'ORG_RES' )
		elapsed_2 = (time.clock() - start_2)
		print "\t\tTime_2: "+str(elapsed_2)

	print "# ---------------------------------------------------------------------------------------------------------------------------------------"
	print "# Run Curvelets: "+DATA_FOLDER
	print "# ---------------------------------------------------------------------------------------------------------------------------------------"
	if runCurv == 1:
		start_1 = time.clock()
		FILE_GFP_BS_CV = FILE_GFP_BS+'_CV'
		runCurv_log = LOCAL_DATASET_PATH_LOG +'/runCurvelets.log'
		TEMP_FILE = open(runCurv_log, 'w')
		TEMP_FILE.write('CurvLog\n')
		TEMP_FILE.close()
		a040_RunCurv.main( FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_PARAMETERS, FILE_GFP_BS, runCurv_log )
		elapsed_1 = (time.clock() - start_1)
		print "\t\tTime_1: "+str(elapsed_1)

	if runCurv_db == 1:
		start_2 = time.clock()
		runCurv_db_log = LOCAL_DATASET_PATH_LOG +'/runCurvelets_db.log'
		TEMP_FILE = open(runCurv_db_log, 'w')
		TEMP_FILE.write('CurvLog\n')
		TEMP_FILE.close()
		a021_Project.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_DATA_DEBUG, FILE_GFP_BS_CV, runCurv_db_log, 'ORG_RES' )
		elapsed_2 = (time.clock() - start_2)
		print "\t\tTime_2: "+str(elapsed_2)

	print "# ---------------------------------------------------------------------------------------------------------------------------------------"
	print "# Run Rescaling: "+DATA_FOLDER
	print "# ---------------------------------------------------------------------------------------------------------------------------------------"
	if runRescale == 1:
		start_1 = time.clock()
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
		elapsed_1 = (time.clock() - start_1)
		print "\t\tTime_1: "+str(elapsed_1)	

	if runRescale_db == 1:
		start_2 = time.clock()
		runRescale_db_log = LOCAL_DATASET_PATH_LOG +'/runRescale_db.log'
		TEMP_FILE = open(runRescale_db_log, 'w')
		TEMP_FILE.write('RescaleLog\n')
		TEMP_FILE.close()

		if haveCy5 == 1:
			a021_Project.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_DATA_DEBUG, FILE_Cy5_BS_RE, runRescale_db_log, 'ORG_RES' )
		if haveTRT == 1:
			a021_Project.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_DATA_DEBUG, FILE_TRI_BS_RE, runRescale_db_log, 'ORG_RES' )
		if haveGFP == 1:
			a021_Project.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_DATA_DEBUG, FILE_GFP_BS_CV_RE, runRescale_db_log, 'ORG_RES' )
		if haveDAP == 1:
			a021_Project.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_DATA_DEBUG, FILE_DAP_BS_RE, runRescale_db_log, 'ORG_RES' )
		elapsed_2 = (time.clock() - start_2)
		print "\t\tTime_2: "+str(elapsed_2)

	print "# ---------------------------------------------------------------------------------------------------------------------------------------"
	print "# Run Rescaling_bit: "+DATA_FOLDER
	print "# ---------------------------------------------------------------------------------------------------------------------------------------"
	if runRescale_bit == 1:
		start_1 = time.clock()
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
		if haveDAP == 1:
			FILE_DAP_BS_RE_bit = FILE_DAP_BS+'_RE_bit'
			a051_RunResc8bit.main( FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_PARAMETERS, FILE_DAP_BS_RE_bit, FILE_DAP_BS, runRescale_bit_log )	
		elapsed_1 = (time.clock() - start_1)
		print "\t\tTime_1: "+str(elapsed_1)

	if runRescale_bit_db == 1:
		start_2 = time.clock()
		runRescale_bit_db_log = LOCAL_DATASET_PATH_LOG +'/runRescale_bit_db.log'
		TEMP_FILE = open(runRescale_bit_db_log, 'w')
		TEMP_FILE.write('Rescale_bitLog\n')
		TEMP_FILE.close()

		if haveCy5 == 1:
			a022_Project8bits.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_DATA_DEBUG, FILE_Cy5_BS_RE_bit, runRescale_bit_db_log, 'ORG_RES' )
		if haveTRT == 1:
			a022_Project8bits.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_DATA_DEBUG, FILE_TRI_BS_RE_bit, runRescale_bit_db_log, 'ORG_RES' )
		if haveGFP == 1:
			a022_Project8bits.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_DATA_DEBUG, FILE_GFP_BS_CV_RE_bit, runRescale_bit_db_log, 'ORG_RES' )
		if haveDAP == 1:
			a022_Project8bits.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_DATA_DEBUG, FILE_DAP_BS_RE_bit, runRescale_bit_db_log, 'ORG_RES' )
		elapsed_2 = (time.clock() - start_2)
		print "\t\tTime_2: "+str(elapsed_2)

	print "# ---------------------------------------------------------------------------------------------------------------------------------------"
	print "# Run Segmentation: "+DATA_FOLDER
	print "# ---------------------------------------------------------------------------------------------------------------------------------------"
	if runSegm == 1:
		start_1 = time.clock()
		runSegm_log = LOCAL_DATASET_PATH_LOG +'/runSegmentation.log'
		TEMP_FILE = open(runSegm_log, 'w')
		TEMP_FILE.write('RunSegmentationLog\n')
		TEMP_FILE.close()

		optionsSegm = LOCAL_DATASET_PATH_SEGM +'/options_segmentation'
		TEMP_FILE = open(optionsSegm, 'w')
		#TEMP_FILE.write('-xTile 700\n')
		#TEMP_FILE.write('-yTile 700\n')
		#TEMP_FILE.write('-zTile 350\n')
		#TEMP_FILE.write('-xTileBor 200\n')
		#TEMP_FILE.write('-yTileBor 200\n')
		#TEMP_FILE.write('-zTileBor 100\n')
		TEMP_FILE.write('-xTile 104\n')
		TEMP_FILE.write('-yTile 104\n')
		TEMP_FILE.write('-zTile 104\n')
		TEMP_FILE.write('-xTileBor 20\n')
		TEMP_FILE.write('-yTileBor 20\n')
		TEMP_FILE.write('-zTileBor 20\n')
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
		TEMP_FILE.close()

		#FILE_Cy5_BS_RE_bit = FILE_Cy5_BS+'_RE_bit'
		a060_RunSegmentation.main( FILE_GFP_BS_CV_RE_bit, FARSIGHT_BIN_EXE, optionsSegm, runSegm_log )
		if( REMOVE_TEMP_SEG == 1 ):
			shutil.rmtree(LOCAL_DATASET_PATH_SEGM_TEMP)
			os.makedirs(LOCAL_DATASET_PATH_SEGM_TEMP)
		elapsed_1 = (time.clock() - start_1)
		print "\t\tTime_1: "+str(elapsed_1)

	if runSegm_db == 1:
		start_2 = time.clock()
		runSegm_db_log = LOCAL_DATASET_PATH_LOG +'/runSegm_db.log'
		TEMP_FILE = open(runSegm_db_log, 'w')
		TEMP_FILE.write('SegmentLog\n')
		TEMP_FILE.close()

		a021_Project.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_SEGM_DEBUG, FILE_GFP_BS_CV_RE_bit+'_soma', runSegm_db_log, 'ORG_RES_BIN' )
		a021_Project.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_SEGM_DEBUG, FILE_GFP_BS_CV_RE_bit+'_label', runSegm_db_log, 'ORG_RES_BIN' )
		a023_ProjectRGB.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_SEGM_DEBUG, FILE_GFP_BS_CV_RE_bit, FILE_GFP_BS_CV_RE_bit+'_soma', '_GFP_SOMA_', runSegm_db_log )
		a023_ProjectRGB.main( LOCAL_DATASET_PATH_LOG, FARSIGHT_BIN_EXE, LOCAL_DATASET_PATH_SEGM_DEBUG, FILE_DAP_BS_RE_bit, FILE_GFP_BS_CV_RE_bit+'_label', '_DAPI_LABEL_', runSegm_db_log )
		elapsed_2 = (time.clock() - start_2)
		print "\t\tTime_2: "+str(elapsed_2)
		


	print "# ---------------------------------------------------------------------------------------------------------------------------------------"
	print "# Run Tracing: "+DATA_FOLDER
	print "# ---------------------------------------------------------------------------------------------------------------------------------------"


	print "# ---------------------------------------------------------------------------------------------------------------------------------------"
	print "# Move Files: "+DATA_FOLDER
	print "# ---------------------------------------------------------------------------------------------------------------------------------------"

	if( MOVE_RESULTS == 1):
		if( MOVE_INTER_STEPS == 0):
			if haveCy5 == 1:
				os.remove(FILE_Cy5+'.nrrd')
				os.remove(FILE_Cy5_BS+'.nrrd')
				os.remove(FILE_Cy5_BS_RE+'.nrrd')
				os.remove(FILE_Cy5_BS_RE_bit+'.nrrd')
			if haveTRT == 1:
				os.remove(FILE_TRI+'.nrrd')
				os.remove(FILE_TRI_BS+'.nrrd')
				os.remove(FILE_TRI_BS_RE+'.nrrd')
				os.remove(FILE_TRI_BS_RE_bit+'.nrrd')
			if haveGFP == 1:
				os.remove(FILE_GFP+'.nrrd')
				os.remove(FILE_GFP_BS+'.nrrd')
				#os.remove(FILE_GFP_BS_CV+'.nrrd')
				os.remove(FILE_GFP_BS_CV_RE+'.nrrd')
				#os.remove(FILE_GFP_BS_CV_RE_bit+'.nrrd')
			if haveDAP == 1:
				os.remove(FILE_DAP+'.nrrd')
				os.remove(FILE_DAP_BS+'.nrrd')
				os.remove(FILE_DAP_BS_RE+'.nrrd')
				os.remove(FILE_DAP_BS_RE_bit+'.nrrd')

			#shutil.rmtree(LOCAL_DATASET_PATH_SEGM_TEMP)

		#shutil.copytree(LOCAL_DATASET_PATH,GLOBAL_DATASET_PATH_RESULTS)

	# ----------------------------------------------

		os.makedirs(GLOBAL_DATASET_PATH_RESULTS+'/Data')
		dirList=os.listdir(LOCAL_DATASET_PATH_DATA)
		for file in dirList:
			f = open(LOCAL_DATASET_PATH_DATA+'/'+file,'r')
			f.close()
			shutil.copy(LOCAL_DATASET_PATH_DATA+'/'+file, GLOBAL_DATASET_PATH_RESULTS+'/Data/'+file)
	

		os.makedirs(GLOBAL_DATASET_PATH_RESULTS+'/Exe')
		dirList=os.listdir(LOCAL_DATASET_PATH_EXE)
		for file in dirList:
			f = open(LOCAL_DATASET_PATH_EXE+'/'+file,'r')
			f.close()
			shutil.copy(LOCAL_DATASET_PATH_EXE+'/'+file, GLOBAL_DATASET_PATH_RESULTS+'/Exe/'+file)

		os.makedirs(GLOBAL_DATASET_PATH_RESULTS+'/Log')
		dirList=os.listdir(LOCAL_DATASET_PATH_LOG)
		for file in dirList:
			f = open(LOCAL_DATASET_PATH_LOG+'/'+file,'r')
			f.close()
			shutil.copy(LOCAL_DATASET_PATH_LOG+'/'+file, GLOBAL_DATASET_PATH_RESULTS+'/Log/'+file)

		os.makedirs(GLOBAL_DATASET_PATH_RESULTS+'/Parameters')
		dirList=os.listdir(LOCAL_DATASET_PATH_PARAMETERS)
		for file in dirList:
			f = open(LOCAL_DATASET_PATH_PARAMETERS+'/'+file,'r')
			f.close()
			shutil.copy(LOCAL_DATASET_PATH_PARAMETERS+'/'+file, GLOBAL_DATASET_PATH_RESULTS+'/Parameters/'+file)

		os.makedirs(GLOBAL_DATASET_PATH_RESULTS+'/Segm')
		dirList=os.listdir(LOCAL_DATASET_PATH_SEGM)
		for file in dirList:
			f = open(LOCAL_DATASET_PATH_SEGM+'/'+file,'r')
			f.close()
			shutil.copy(LOCAL_DATASET_PATH_SEGM+'/'+file, GLOBAL_DATASET_PATH_RESULTS+'/Segm/'+file)

		os.makedirs(GLOBAL_DATASET_PATH_RESULTS+'/TracesAndSomas')
		dirList=os.listdir(LOCAL_DATASET_PATH_TRACE_SOMAS)
		for file in dirList:
			f = open(LOCAL_DATASET_PATH_TRACE_SOMAS+'/'+file,'r')
			f.close()
			shutil.copy(LOCAL_DATASET_PATH_TRACE_SOMAS+'/'+file, GLOBAL_DATASET_PATH_RESULTS+'/TracesAndSomas/'+file)
	# ----------------------------------------------


		#TEMP2 = 'cp -r '+LOCAL_DATASET_PATH_DATA+' '+GLOBAL_DATASET_PATH_RESULTS+'/'
		#print TEMP2
		#TEMP3 = subprocess.Popen(TEMP2, shell=True)
		##subprocess.call(['cp -r '])
		#TEMP3.communicate()
		#TEMP2 = 'cp -r '+LOCAL_DATASET_PATH_EXE+' '+GLOBAL_DATASET_PATH_RESULTS+'/'
		#print TEMP2
		#TEMP4 = subprocess.Popen(TEMP2, shell=True)
		#TEMP4.communicate()
		#TEMP2 = 'cp -r '+LOCAL_DATASET_PATH_LOG+' '+GLOBAL_DATASET_PATH_RESULTS+'/'
		#print TEMP2
		#TEMP5 = subprocess.Popen(TEMP2, shell=True)
		#TEMP5.communicate()
		#TEMP2 = 'cp -r '+LOCAL_DATASET_PATH_PARAMETERS+' '+GLOBAL_DATASET_PATH_RESULTS+'/'
		#print TEMP2
		#TEMP6 = subprocess.Popen(TEMP2, shell=True)
		#TEMP6.communicate()
		#TEMP2 = 'cp -r '+LOCAL_DATASET_PATH_SEGM+' '+GLOBAL_DATASET_PATH_RESULTS+'/'
		#print TEMP2
		#TEMP7 = subprocess.Popen(TEMP2, shell=True)
		#TEMP7.communicate()
		#TEMP2 = 'cp -r '+LOCAL_DATASET_PATH_TRACE_SOMAS+' '+GLOBAL_DATASET_PATH_RESULTS+'/'
		#print TEMP2
		#TEMP8 = subprocess.Popen(TEMP2, shell=True)
		#TEMP8.communicate()

		#TEMP2 = 'rm -rf '+LOCAL_DATASET_PATH_DEBUG
		#TEMP3 = subprocess.Popen(TEMP2, shell=True)
		#TEMP3.communicate()

		#shutil.copytree(LOCAL_DATASET_PATH_DATA,GLOBAL_DATASET_PATH_RESULTS+'/Data') #// Genera errores
		#shutil.copytree(LOCAL_DATASET_PATH_EXE,GLOBAL_DATASET_PATH_RESULTS+'/Exe')
		#shutil.copytree(LOCAL_DATASET_PATH_LOG,GLOBAL_DATASET_PATH_RESULTS+'/Log')
		#shutil.copytree(LOCAL_DATASET_PATH_PARAMETERS,GLOBAL_DATASET_PATH_RESULTS+'/Parameters')
		#shutil.copytree(LOCAL_DATASET_PATH_SEGM,GLOBAL_DATASET_PATH_RESULTS+'/Segm')
		#shutil.copytree(LOCAL_DATASET_PATH_TRACE_SOMAS,GLOBAL_DATASET_PATH_RESULTS+'/TracesAndSomas')

	if( REMOVE_MONTAGES == 1):
		shutil.rmtree(LOCAL_DATASET_PATH_DATA)
		
	print "# ---------------------------------------------------------------------------------------------------------------------------------------"
	print "# Move internally: "+DATA_FOLDER
	print "# ---------------------------------------------------------------------------------------------------------------------------------------"
	if( MOVE_LOCALLY == 1):
		shutil.rmtree(LOCAL_DEB_DATASET_PATH)
		shutil.copytree(LOCAL_DATASET_PATH,LOCAL_DEB_DATASET_PATH) #// Genera errores
		shutil.rmtree(LOCAL_DATASET_PATH)

	#print "# ---------------------------------------------------------------------------------------------------------------------------------------"
	#print "# Test if the results are the same as the previous ones: "+DATA_FOLDER
	#print "# ---------------------------------------------------------------------------------------------------------------------------------------"
