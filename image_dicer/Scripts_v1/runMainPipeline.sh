#!/bin/bash
# This file is the batch script that will run the PIPELINE of DARPA datasets

# ##############################################################################################################################
# # Make Farsight
# ##############################################################################################################################
cd $FARSIGHT_BIN
make -j80
cd $ACTUAL_DIRECTORY


# ##############################################################################################################################
# # Move Images
# ##############################################################################################################################
/usr/bin/time $ACTUAL_DIRECTORY/runMoveImages.sh > $LOCAL_DATASET_PATH_LOG/runMoveImages.log 2>&1

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
/usr/bin/time $ACTUAL_DIRECTORY/runBackgroundsubstraction.sh > $LOCAL_DATASET_PATH_LOG/runBackgroundsubstraction.log 2>&1


# ##############################################################################################################################
# # Curvelets GFP channel
# ##############################################################################################################################
cp $FARSIGHT_BIN_EXE/curvelets $LOCAL_DATASET_PATH_EXE
/usr/bin/time $ACTUAL_DIRECTORY/runCurvelets.sh > $LOCAL_DATASET_PATH_LOG/runCurvelets.log 2>&1


# ##############################################################################################################################
# # Segmentation, features computation, soma extraction
# ##############################################################################################################################
/usr/bin/time $ACTUAL_DIRECTORY/runSegmentation.sh > $LOCAL_DATASET_PATH_LOG/runSegmentation.log #2>&1
 
 
##############################################################################################################################
# # Tracing
##############################################################################################################################
/usr/bin/time $ACTUAL_DIRECTORY/runTracing.sh > $LOCAL_DATASET_PATH_LOG/runTracing.log 2>&1

##############################################################################################################################
# # Copy results to far-01
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
