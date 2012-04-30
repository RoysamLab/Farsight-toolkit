#!/bin/bash
echo $FARSIGHT_BIN_EXE/darpa_tracer_w_seg_v2 $DAPI_LOCAL\_BS.nrrd $GFP_LOCAL\_BS_CV.mhd $Cy5_LOCAL\_BS.nrrd $LOCAL_PARAMETERS_PATH/Seg_Params.ini $LOCAL_PARAMETERS_PATH/ProjectDefinition.xml 0 $SMALLIMAGE 1 $LOCAL_PARAMETERS_PATH/options_mnt

# Segmentation part 1
$FARSIGHT_BIN_EXE/darpa_tracer_w_seg_v2 $DAPI_LOCAL\_BS.nrrd $GFP_LOCAL\_BS_CV.mhd $Cy5_LOCAL\_BS.nrrd $LOCAL_PARAMETERS_PATH/Seg_Params.ini $LOCAL_PARAMETERS_PATH/ProjectDefinition.xml 0 $SMALLIMAGE 2 $LOCAL_PARAMETERS_PATH/options_mnt
$FARSIGHT_BIN_EXE/darpa_tracer_w_seg_v2 $DAPI_LOCAL\_BS.nrrd $GFP_LOCAL\_BS_CV.mhd $Cy5_LOCAL\_BS.nrrd $LOCAL_PARAMETERS_PATH/Seg_Params.ini $LOCAL_PARAMETERS_PATH/ProjectDefinition.xml 0 $SMALLIMAGE 3 $LOCAL_PARAMETERS_PATH/options_mnt

# Segmentation Stiching part 2
# echo segmentation part 2
$FARSIGHT_BIN_EXE/darpa_tracer_w_seg_v2 $DAPI_LOCAL\_BS.nrrd $GFP_LOCAL\_BS_CV.mhd $Cy5_LOCAL\_BS.nrrd $LOCAL_PARAMETERS_PATH/Seg_Params.ini $LOCAL_PARAMETERS_PATH/ProjectDefinition.xml 4 $SMALLIMAGE 1 $LOCAL_PARAMETERS_PATH/options_mnt



#/data/research/Fiji.app/fiji-linux64 --headless -macro $LOCAL_PARAMETERS_PATH/fijiMacro_maxProj.ijm $GFP_LOCAL\_BS_CV_soma_montage -batch


for f in $GFP_LOCAL\_BS_CV.mhd
do
	export GFP_LOCAL_BS_CV_EXE=$f
done
export GFP_LOCAL_BS_CV=${GFP_LOCAL_BS_CV_EXE%\.*}
/data/research/Fiji.app/fiji-linux64 --headless -macro $LOCAL_PARAMETERS_PATH/fijiMacro_maxProj.ijm $GFP_LOCAL_BS_CV -batch

for f in $GFP_LOCAL\_BS_CV_soma_montage.mhd
do
	export GFP_LOCAL_BS_CV_soma_montage_EXE=$f
done
export GFP_LOCAL_BS_CV_soma_montage=${GFP_LOCAL_BS_CV_soma_montage_EXE%\.*}
/data/research/Fiji.app/fiji-linux64 --headless -macro $LOCAL_PARAMETERS_PATH/fijiMacro_maxProj.ijm $GFP_LOCAL_BS_CV_soma_montage -batch

# LOCAL_GFP_MHD_EXT=$LOCAL_DATASET_PATH/*GFPdsu.mhd
# for f in $LOCAL_GFP_MHD_EXT
# do
# 	export GFP_LOCAL_EXE=$f
# done
# export GFP_LOCAL=${GFP_LOCAL_EXE%\.*}

# onlyTrace $SMALLIMAGE segSteps
