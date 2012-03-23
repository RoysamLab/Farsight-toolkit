#!/bin/bash
# This file is the batch script that will run the PIPELINE of DARPA datasets

##############################################################################################################################
# PARAMETERS
##############################################################################################################################
ACTUAL_DIRECTORY=$(pwd)

# The directory where the data is stored in the local machine (FAR-04, FAR-05)
# export DATA_FOLDER=0117\_Tile
export FARSIGHT_BIN=/data/nicolas/farsight_updated/bin
export FARSIGHT_BIN_EXE=/data/nicolas/farsight_updated/bin/exe

export LOCAL_DATASET_PATH=/data/nicolas/data/$DATA_FOLDER					# ---> This directory has to exist
export LOCAL_PARAMETERS_PATH=/data/nicolas/data/$DATA_FOLDER/Parameters		# ---> This directory has to exist
export LOCAL_DATASET_PATH_EXE=/data/nicolas/data/$DATA_FOLDER/Exe
export LOCAL_DATASET_PATH_LOG=/data/nicolas/data/$DATA_FOLDER/Log
export LOCAL_DATASET_PATH_TRACE_SOMAS=/data/nicolas/data/$DATA_FOLDER/TracesAndSomas

export GLOBAL_DATASET_PATH=/FSdata/data/$DATA_FOLDER						# ---> This directory has to exist
export GLOBAL_DATASET_PATH_RESULTS=/FSdata/data/$DATA_FOLDER\_RESULTS_FAR04


export REMOVE_MONTAGES=0 	# This flag is set in case we want the montages to be removed after the process is done, especially when running many montages in serial we want to make sure not to run out of memory
export MOVE_RESULTS=3		# If 1 the results will be moved
			# if 0 the results will be keep
			# if 2 the results will be copied (keep and copy to FSDATA)
			# if 3 move everysingle file, exept the folder, which are copied

if [ ! -d $LOCAL_DATASET_PATH_EXE ]; then
	mkdir $LOCAL_DATASET_PATH_EXE
else
	echo L1n: The content of the path $LOCAL_DATASET_PATH_EXE will be erased
	cd $LOCAL_DATASET_PATH_EXE
# 	rm -rf $LOCAL_DATASET_PATH_EXE/*
	cd $ACTUAL_DIRECTORY
fi

if [ ! -d $LOCAL_DATASET_PATH_TRACE_SOMAS ]; then
	mkdir $LOCAL_DATASET_PATH_TRACE_SOMAS
else
	echo L1n: The content of the path $LOCAL_DATASET_PATH_TRACE_SOMAS will be erased
	cd $LOCAL_DATASET_PATH_TRACE_SOMAS
# 	rm -rf $LOCAL_DATASET_PATH_EXE/*
	cd $ACTUAL_DIRECTORY
fi

# if [ ! -d $DATASET_PATH_RESULTS ]; then
# 	mkdir $DATASET_PATH_RESULTS
# else
# 	echo L1n: The content of the path $DATASET_PATH_RESULTS will be erased
# 	cd $DATASET_PATH_RESULTS
# # 	rm -rf $DATASET_PATH_RESULTS/*
# 	cd $ACTUAL_DIRECTORY
# fi

if [ ! -d $LOCAL_DATASET_PATH_LOG ]; then
	mkdir $LOCAL_DATASET_PATH_LOG
else
	echo L1n: The content of the path $LOCAL_DATASET_PATH_LOG will be erased
	cd $LOCAL_DATASET_PATH_LOG
# 	rm -rf $DATASET_PATH_RESULTS/*
	cd $ACTUAL_DIRECTORY
fi


# automatically export EVERYTHING
# set -a
