#!/bin/bash
# This file is the batch script that will run the PIPELINE of DARPA datasets

##############################################################################################################################
# PARAMETERS
##############################################################################################################################
ACTUAL_DIRECTORY=$(pwd)

# The directory where the data is stored in the local machine (FAR-04, FAR-05)
export DATA_FOLDER=0117_Tile
export LOCAL_DATASET_PATH=/data/nicolas/data/$DATA_FOLDER			# ---> This directory has to exist
export LOCAL_PARAMETERS_PATH=/data/nicolas/data/$DATA_FOLDER/Parameters		# ---> This directory has to exist
export LOCAL_DATASET_PATH_EXE=/data/nicolas/data/$DATA_FOLDER/Exe
export LOCAL_DATASET_PATH_LOG=/data/nicolas/data/$DATA_FOLDER/Log
export LOCAL_DATASET_PATH_TRACE_SOMAS=/data/nicolas/data/$DATA_FOLDER/TracesAndSomas

export GLOBAL_DATASET_PATH=/FSdata/data/$DATA_FOLDER				# ---> This directory has to exist
export GLOBAL_DATASET_PATH_RESULTS=/FSdata/data/$DATA_FOLDER\_RESULTS


export FARSIGHT_BIN=/data/nicolas/farsigt3_from05/bin
export FARSIGHT_BIN_EXE=/data/nicolas/farsigt3_from05/bin/exe

export REMOVE_MONTAGES=0 	# This flag is set in case we want the montages to be removed after the process is done, especially when running many montages in serial we want to make sure not to run out of memory
export MOVE_RESULTS=3		# If 1 the results will be moved
			# if 0 the results will be keep
			# if 2 the results will be copied (keep and copy to FSDATA)
# if 3 move everysingle file, exept the folder, which are copied



# Curvelets PATH
export CURVELETS_PATH=/data/nicolas/farsigt3_from05/bin/exe

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


# echo $n

# ##############################################################################################################################
# # Make Farsight
# ##############################################################################################################################
cd $FARSIGHT_BIN
make -j80
cd $ACTUAL_DIRECTORY


# ##############################################################################################################################
# # Move Images
# ##############################################################################################################################
/usr/bin/time $ACTUAL_DIRECTORY/moveImages.sh -> $LOCAL_DATASET_PATH_LOG/moveImages.log

LOCAL_DAPI_MHD_EXT=$LOCAL_DATASET_PATH/*DAPIdsu.mhd
for f in $LOCAL_DAPI_MHD_EXT
do
	export DAPI_LOCAL_EXE=$f
done
export DAPI_LOCAL=${DAPI_LOCAL_EXE%\.*}

LOCAL_Cy5_MHD_EXT=$LOCAL_DATASET_PATH/*Cy5dsu.mhd
for f in $LOCAL_Cy5_MHD_EXT
do
	export Cy5_LOCAL_EXE=$f
done
export Cy5_LOCAL=${Cy5_LOCAL_EXE%\.*}

LOCAL_TRITC_MHD_EXT=$LOCAL_DATASET_PATH/*TRITCdsu.mhd
for f in $LOCAL_TRITC_MHD_EXT
do
	export TRITC_LOCAL_EXE=$f
done
export TRITC_LOCAL=${TRITC_LOCAL_EXE%\.*}

LOCAL_GFP_MHD_EXT=$LOCAL_DATASET_PATH/*GFPdsu.mhd
for f in $LOCAL_GFP_MHD_EXT
do
	export GFP_LOCAL_EXE=$f
done
export GFP_LOCAL=${GFP_LOCAL_EXE%\.*}



# ##############################################################################################################################
# # Run Background Substraction
# ##############################################################################################################################
# echo $ACTUAL_DIRECTORY
/usr/bin/time $ACTUAL_DIRECTORY/runBackgroundsubstraction.sh > $LOCAL_DATASET_PATH_LOG/runBackgroundsubstraction.log


# ##############################################################################################################################
# # # Curvelets GFP channel
# ##############################################################################################################################
cp $CURVELETS_PATH/curvelets $LOCAL_DATASET_PATH_EXE
/usr/bin/time $ACTUAL_DIRECTORY/runCurvelets.sh > $LOCAL_DATASET_PATH_LOG/runCurvelets.log


# # ##############################################################################################################################
# # # Segmentation, features computation, soma extraction
# # ##############################################################################################################################
# cd $FARSIGHT_BIN_EXE
# echo $DAPI_LOCAL\_BS.nrrd
# echo $GFP_LOCAL\_BS._CV.mhd
# echo $Cy5_LOCAL\_BS.nrrd
# echo $LOCAL_PARAMETERS_PATH/ProjectDefinition.xml
./$FARSIGHT_BIN_EXE/darpa_tracer_w_seg $DAPI_LOCAL\_BS.nrrd $GFP_LOCAL\_BS_CV.mhd $Cy5_LOCAL\_BS.nrrd $LOCAL_PARAMETERS_PATH/Seg_Params.ini $LOCAL_PARAMETERS_PATH/ProjectDefinition.xml 0 > $LOCAL_DATASET_PATH_LOG/runSegmentation.log
 
 
##############################################################################################################################
# Tracing
##############################################################################################################################
# cd $FARSIGHT_BIN_EXE
$FARSIGHT_BIN_EXE/darpa_tracer_w_seg $DAPI_LOCAL\_BS.nrrd $GFP_LOCAL\_BS_CV.mhd $Cy5_LOCAL\_BS.nrrd $LOCAL_PARAMETERS_PATH/Seg_Params.ini $LOCAL_PARAMETERS_PATH/ProjectDefinition.xml 1 > $LOCAL_DATASET_PATH_LOG/runTracing.log


##############################################################################################################################
# Copy results to far-01
##############################################################################################################################
if [ ! -d $GLOBAL_DATASET_PATH_RESULTS ]; then
	echo L1n: Result Path does not exist, it will be created
	mkdir $GLOBAL_DATASET_PATH_RESULTS
# 	cd $GLOBAL_DATASET_PATH_RESULTS
# else
# 	echo L1n: The content of the path $DATASET_PATH_RESULTS will be erased
# 	cd $DATASET_PATH_RESULTS
# 	rm -rf $DATASET_PATH_RESULTS *
fi

# export MOVE_RESULTS=2		# If 1 the results will be moved
			# if 0 the results will be keep
			# if 2 the results will be copied (keep and copy to FSDATA)


# echo $GLOBAL_DATASET_PATH_RESULTS
if [ $MOVE_RESULTS -eq "2" ]; then
	cd $LOCAL_DATASET_PATH
# 	echo si
# 	cp *label.* 	$GLOBAL_DATASET_PATH_RESULTS
# 	cp *soma_* 	$GLOBAL_DATASET_PATH_RESULTS
# 	cp OnlySWC.xml	$GLOBAL_DATASET_PATH_RESULTS
# 	cp *.swc 	$GLOBAL_DATASET_PATH_RESULTS
# 	cp * 		$GLOBAL_DATASET_PATH_RESULTS # Copy all the files and folders
	cp -r * 	$GLOBAL_DATASET_PATH_RESULTS # Copy all the files
	cd $ACTUAL_DIRECTORY
# else
# 	echo no
fi
if [ $MOVE_RESULTS -eq "1" ]; then
	cd $LOCAL_DATASET_PATH
# 	echo si
# 	mv *label.* 	$GLOBAL_DATASET_PATH_RESULTS
# 	mv *soma_* 	$GLOBAL_DATASET_PATH_RESULTS
# 	mv OnlySWC.xml	$GLOBAL_DATASET_PATH_RESULTS
# 	mv *.swc 	$GLOBAL_DATASET_PATH_RESULTS
	mv -r * 	$GLOBAL_DATASET_PATH_RESULTS # Copy all the files
	cd $ACTUAL_DIRECTORY
# else
# 	echo no
fi
if [ $MOVE_RESULTS -eq "3" ]; then
	cd $LOCAL_DATASET_PATH
# 	echo si
# 	mv *label.* 	$GLOBAL_DATASET_PATH_RESULTS
# 	mv *soma_* 	$GLOBAL_DATASET_PATH_RESULTS
# 	mv OnlySWC.xml	$GLOBAL_DATASET_PATH_RESULTS
# 	mv *.swc 	$GLOBAL_DATASET_PATH_RESULTS
	mv *.* 		$GLOBAL_DATASET_PATH_RESULTS # Copy all the files
	cp -r * 	$GLOBAL_DATASET_PATH_RESULTS # Copy all the files
	cd $ACTUAL_DIRECTORY
# else
# 	echo no
fi
