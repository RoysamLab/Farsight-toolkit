#!/bin/bash

echo About to run curvelets
$LOCAL_DATASET_PATH_EXE/curvelets $GFP_LOCAL'_BS.nrrd' $LOCAL_PARAMETERS_PATH/options_curvelets
# $LOCAL_DATASET_PATH_EXE/curvelets $GFP_LOCAL'_BS2.nrrd' $LOCAL_PARAMETERS_PATH/options_curvelets

for f in $GFP_LOCAL\_BS_CV.mhd
do
	export GFP_LOCAL_BS_CV_EXE=$f
done
export GFP_LOCAL_BS_CV=${GFP_LOCAL_BS_CV_EXE%\.*}
/data/research/Fiji.app/fiji-linux64 --headless -macro $LOCAL_PARAMETERS_PATH/fijiMacro_maxProj.ijm $GFP_LOCAL_BS_CV -batch