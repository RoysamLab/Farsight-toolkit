#!/bin/bash

# This file is the batch script that will run the PIPELINE of DARPA datasets

# dont forget this
# chmod u+x ../../farsigt3_from05/src/image_dicer/*_v2.sh
# . Source to export back into the parent


# export DATA_FOLDER=0113
# . ./runCreateFolders_v2.sh
# . ./runMainPipeline_v2.sh > $LOCAL_DATASET_PATH_LOG/runMainPipeline.log
# export DATA_FOLDER=0120
# . ./runCreateFolders_v2.sh
# . ./runMainPipeline_v2.sh > $LOCAL_DATASET_PATH_LOG/runMainPipeline.log
# export DATA_FOLDER=0117
# . ./runCreateFolders_v2.sh
# . ./runMainPipeline_v2.sh > $LOCAL_DATASET_PATH_LOG/runMainPipeline.log
# export DATA_FOLDER=0123
# . ./runCreateFolders_v2.sh
# . ./runMainPipeline_v2.sh > $LOCAL_DATASET_PATH_LOG/runMainPipeline.log
# export DATA_FOLDER=0103
# . ./runCreateFolders_v2.sh
# . ./runMainPipeline_v2.sh > $LOCAL_DATASET_PATH_LOG/runMainPipeline.log
#export DATA_FOLDER=0102
#. ./runCreateFolders_v2.sh
#. ./runMainPipeline_v2.sh > $LOCAL_DATASET_PATH_LOG/runMainPipeline.log
#export DATA_FOLDER=0131
#. ./runCreateFolders_v2.sh
#. ./runMainPipeline_v2.sh > $LOCAL_DATASET_PATH_LOG/runMainPipeline.log # 2>&1
# export DATA_FOLDER=0128
# . ./runCreateFolders_v2.sh
# . ./runMainPipeline_v2.sh > $LOCAL_DATASET_PATH_LOG/runMainPipeline.log



# export DATA_FOLDER=0102
# . ./runCreateFolders_v2.sh
# . ./runMainPipeline_v2.sh > $LOCAL_DATASET_PATH_LOG/runMainPipeline.log # 2>&1

# export DATA_FOLDER=0103
# . ./runCreateFolders_v2.sh
# . ./runMainPipeline_v2.sh > $LOCAL_DATASET_PATH_LOG/runMainPipeline.log
# export DATA_FOLDER=0113
# . ./runCreateFolders_v2.sh
# . ./runMainPipeline_v2.sh > $LOCAL_DATASET_PATH_LOG/runMainPipeline.log
 export DATA_FOLDER=0117 #_Tile_mhd
 . ./runCreateFolders_v2.sh
 . ./runMainPipeline_v2_debug.sh > $LOCAL_DATASET_PATH_LOG/runMainPipeline.log
# export DATA_FOLDER=0120
# . ./runCreateFolders_v2.sh
# . ./runMainPipeline_v2.sh > $LOCAL_DATASET_PATH_LOG/runMainPipeline.log
# export DATA_FOLDER=0123
# . ./runCreateFolders_v2.sh
# . ./runMainPipeline_v2.sh > $LOCAL_DATASET_PATH_LOG/runMainPipeline.log
# export DATA_FOLDER=0128
# . ./runCreateFolders_v2.sh
# . ./runMainPipeline_v2.sh > $LOCAL_DATASET_PATH_LOG/runMainPipeline.log
# export DATA_FOLDER=0131
# . ./runCreateFolders_v2.sh
# . ./runMainPipeline_v2.sh > $LOCAL_DATASET_PATH_LOG/runMainPipeline.log

# # ./pipelineSegAndTrace_0103_v2.sh
# # ./pipelineSegAndTrace_0113_v2.sh #/far04
# # ./pipelineSegAndTrace_0117_v2.sh #/far05
# # ./pipelineSegAndTrace_0120_v2.sh #/far05
# # ./pipelineSegAndTrace_0123_v2.sh #/far05
# ./pipelineSegAndTrace_0128_v2.sh #/far05
# # ./pipelineSegAndTrace_0131_v2.sh #/far04
# ./pipelineSegAndTrace_1118_v2.sh
# ./pipelineSegAndTrace_1122_v2.sh
# ./pipelineSegAndTrace_1230_v2.sh
