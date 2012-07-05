#!/bin/bash
# function testing(){
if [ -f $LOCAL_DATASET_PATH/*DAPIdsu.raw ]; then
	echo L1n: File exists now we are going to compare it DAPI
	cmp -s $GLOBAL_DATASET_PATH/*DAPIdsu.raw $LOCAL_DATASET_PATH/*DAPIdsu.raw
	if [ $? -eq 1 ]; then
		echo L1n: Copying DAPI montage it was different
		cp -f $GLOBAL_DATASET_PATH/*DAPIdsu.raw $LOCAL_DATASET_PATH
	fi
else
	echo L1n: File does not exists now we are going to copy it DAPI
	cp -f $GLOBAL_DATASET_PATH/*DAPIdsu.raw $LOCAL_DATASET_PATH
fi

if [ -f $LOCAL_DATASET_PATH/*DAPIdsu.mhd ]; then
	echo L1n: File exists now we are going to compare it DAPI
	cmp -s $GLOBAL_DATASET_PATH/*DAPIdsu.mhd $LOCAL_DATASET_PATH/*DAPIdsu.mhd
	if [ $? -eq 1 ]; then
		echo L1n: Copying DAPI montage it was different
		cp -f $GLOBAL_DATASET_PATH/*DAPIdsu.mhd $LOCAL_DATASET_PATH
	fi
else
	echo L1n: File does not exists now we are going to copy it DAPI
	cp -f $GLOBAL_DATASET_PATH/*DAPIdsu.mhd $LOCAL_DATASET_PATH
fi

# if [ -f $LOCAL_DATASET_PATH/*Cy5dsu.raw ]; then
# 	echo L1n: File exists now we are going to compare it Cy5
# 	cmp -s $GLOBAL_DATASET_PATH/*Cy5dsu.raw $LOCAL_DATASET_PATH/*Cy5dsu.raw
# 	if [ $? -eq 1 ]; then
# 		echo L1n: Copying Cy5 montage it was different
# 		cp -f $GLOBAL_DATASET_PATH/*Cy5dsu.raw $LOCAL_DATASET_PATH
# 	fi
# else
# 	echo L1n: File does not exists now we are going to copy it Cy5
# 	cp -f $GLOBAL_DATASET_PATH/*Cy5dsu.raw $LOCAL_DATASET_PATH
# fi
# 
# if [ -f $LOCAL_DATASET_PATH/*Cy5dsu.mhd ]; then
# 	echo L1n: File exists now we are going to compare it Cy5
# 	cmp -s $GLOBAL_DATASET_PATH/*Cy5dsu.mhd $LOCAL_DATASET_PATH/*Cy5dsu.mhd
# 	if [ $? -eq 1 ]; then
# 		echo L1n: Copying Cy5 montage it was different
# 		cp -f $GLOBAL_DATASET_PATH/*Cy5dsu.mhd $LOCAL_DATASET_PATH
# 	fi
# else
# 	echo L1n: File does not exists now we are going to copy it Cy5
# 	cp -f $GLOBAL_DATASET_PATH/*Cy5dsu.mhd $LOCAL_DATASET_PATH
# fi

# if [ -f $LOCAL_DATASET_PATH/*TRITCdsu.raw ]; then
# 	echo L1n: File exists now we are going to compare it TRITC
# 	cmp -s $GLOBAL_DATASET_PATH/*TRITCdsu.raw $LOCAL_DATASET_PATH/*TRITCdsu.raw
# 	if [ $? -eq 1 ]; then
# 		echo L1n: Copying TRITC montage it was different
# 		cp -f $GLOBAL_DATASET_PATH/*TRITCdsu.raw $LOCAL_DATASET_PATH
# 	fi
# else
# 	echo L1n: File does not exists now we are going to copy it TRITC
# 	cp -f $GLOBAL_DATASET_PATH/*TRITCdsu.raw $LOCAL_DATASET_PATH
# fi
# 
# if [ -f $LOCAL_DATASET_PATH/*TRITCdsu.mhd ]; then
# 	echo L1n: File exists now we are going to compare it TRITC
# 	cmp -s $GLOBAL_DATASET_PATH/*TRITCdsu.mhd $LOCAL_DATASET_PATH/*TRITCdsu.mhd
# 	if [ $? -eq 1 ]; then
# 		echo L1n: Copying TRITC montage it was different
# 		cp -f $GLOBAL_DATASET_PATH/*TRITCdsu.mhd $LOCAL_DATASET_PATH
# 	fi
# else
# 	echo L1n: File does not exists now we are going to copy it TRITC
# 	cp -f $GLOBAL_DATASET_PATH/*TRITCdsu.mhd $LOCAL_DATASET_PATH
# fi

if [ -f $LOCAL_DATASET_PATH/*GFPdsu.raw ]; then
	echo L1n: File exists now we are going to compare it GFP
	cmp -s $GLOBAL_DATASET_PATH/*GFPdsu.raw $LOCAL_DATASET_PATH/*GFPdsu.raw
	if [ $? -eq 1 ]; then
		echo L1n: Copying GFP montage it was different
		cp -f $GLOBAL_DATASET_PATH/*GFPdsu.raw $LOCAL_DATASET_PATH
	fi
else
	echo L1n: File does not exists now we are going to copy it GFP
	cp -f $GLOBAL_DATASET_PATH/*GFPdsu.raw $LOCAL_DATASET_PATH
fi

if [ -f $LOCAL_DATASET_PATH/*GFPdsu.mhd ]; then
	echo L1n: File exists now we are going to compare it GFP
	cmp -s $GLOBAL_DATASET_PATH/*GFPdsu.mhd $LOCAL_DATASET_PATH/*GFPdsu.mhd
	if [ $? -eq 1 ]; then
		echo L1n: Copying GFP montage it was different
		cp -f $GLOBAL_DATASET_PATH/*GFPdsu.mhd $LOCAL_DATASET_PATH
	fi
else
	echo L1n: File does not exists now we are going to copy it GFP
	cp -f $GLOBAL_DATASET_PATH/*GFPdsu.mhd $LOCAL_DATASET_PATH
fi

# wait
# }

# /usr/bin/time testing



