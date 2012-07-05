#!/bin/bash

# This file is the batch script that will run the PIPELINE of DARPA datasets

# dont forget this
# chmod u+x ../../farsigt3_from05/src/image_dicer/*_v2.sh
# . Source to export back into the parent

export REMOVE_MONTAGES=0 	# This flag is set in case we want the montages to be removed after the process is done, especially when running many montages in serial we want to make sure not to run out of memory
export MOVE_RESULTS=2		# If 1 the results will be moved
				# if 0 the results will be keep
				# if 2 the results will be copied (keep and copy to FSDATA)
				# if 3 move everysingle file, exept the folder, which are copied
export SMALLIMAGE=1		#

export runMove=1
export runBack=1
export runCurv=1
export runSegm=1
export runTrac=1

# if $SMALLIMAGE
# then


export DATA_FOLDER=0131_crop
. ./runCreateFolders_v2.sh
. ./runMainPipeline_v2.sh > $LOCAL_DATASET_PATH_LOG/runMainPipeline.log

# export DATA_FOLDER=0405
# . ./runCreateFolders_v2.sh
# . ./runMainPipeline_v2.sh > $LOCAL_DATASET_PATH_LOG/runMainPipeline.log

# export DATA_FOLDER=0323 #_crop
# . ./runCreateFolders_v2.sh
# . ./runMainPipeline_v2.sh > $LOCAL_DATASET_PATH_LOG/runMainPipeline.log

# export DATA_FOLDER=0410 #_crop
# . ./runCreateFolders_v2.sh
# . ./runMainPipeline_v2.sh > $LOCAL_DATASET_PATH_LOG/runMainPipeline.log

# export DATA_FOLDER=0405 #_crop
# . ./runCreateFolders_v2.sh
# . ./runMainPipeline_v2.sh > $LOCAL_DATASET_PATH_LOG/runMainPipeline.log
