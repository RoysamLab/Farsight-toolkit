#!/bin/bash

# echo $GFP_LOCAL'_BS.nrrd'

if [ -f $GFP_LOCAL'_BS.nrrd' ]; then
	echo L1n: GFP BS already exist
else
	echo L1n: GFP BS does not exist
	echo Doing backsubstraction on GFP
	/data/research/Fiji.app/fiji-linux64 --headless -macro $LOCAL_PARAMETERS_PATH/fijiMacro.ijm $GFP_LOCAL -batch	
fi

if [ -f $DAPI_LOCAL'_BS.nrrd' ]; then
	echo L1n: DAPI BS already exist
else
	echo L1n: DAPI BS does not exist
	echo Doing backsubstraction on DAPI
	/data/research/Fiji.app/fiji-linux64 --headless -macro $LOCAL_PARAMETERS_PATH/fijiMacro.ijm $DAPI_LOCAL -batch	
fi

# if [ -f $Cy5_LOCAL'_BS.nrrd' ]; then
# 	echo L1n: Cy5 BS already exist
# else
# 	echo L1n: Cy5 BS does not exist
# 	echo Doing backsubstraction on Cy5
# 	/data/research/Fiji.app/fiji-linux64 --headless -macro $LOCAL_PARAMETERS_PATH/fijiMacro.ijm $Cy5_LOCAL -batch	
# fi

# if [ -f $TRITC_LOCAL'_BS.nrrd' ]; then
# 	echo L1n: TRITC BS already exist
# else
# 	echo L1n: TRITC BS does not exist
# 	echo Doing backsubstraction on TRITC
# 	/data/research/Fiji.app/fiji-linux64 --headless -macro $LOCAL_PARAMETERS_PATH/fijiMacro.ijm $TRITC_LOCAL -batch	
# fi

#wait not a good idea to run this in parallel
