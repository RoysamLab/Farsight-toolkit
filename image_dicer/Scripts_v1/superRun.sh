#!/bin/bash

# This file is the batch script that will run the PIPELINE of DARPA datasets

# dont forget this
# chmod u+x ../../farsigt3_from05/src/image_dicer/*.sh
# . Source to export back into the parent

# export DATA_FOLDER=0102
# . ./runCreateFolders.sh
# . ./runMainPipeline.sh > $LOCAL_DATASET_PATH_LOG/runMainPipeline.log # 2>&1
# export DATA_FOLDER=0103
# . ./runCreateFolders.sh
# . ./runMainPipeline.sh > $LOCAL_DATASET_PATH_LOG/runMainPipeline.log
export DATA_FOLDER=0113
. ./runCreateFolders.sh
. ./runMainPipeline.sh > $LOCAL_DATASET_PATH_LOG/runMainPipeline.log
# export DATA_FOLDER=0117
# . ./runCreateFolders.sh
# . ./runMainPipeline.sh > $LOCAL_DATASET_PATH_LOG/runMainPipeline.log
export DATA_FOLDER=0120
. ./runCreateFolders.sh
. ./runMainPipeline.sh > $LOCAL_DATASET_PATH_LOG/runMainPipeline.log
# export DATA_FOLDER=0123
# . ./runCreateFolders.sh
# . ./runMainPipeline.sh > $LOCAL_DATASET_PATH_LOG/runMainPipeline.log
# export DATA_FOLDER=0128
# . ./runCreateFolders.sh
# . ./runMainPipeline.sh > $LOCAL_DATASET_PATH_LOG/runMainPipeline.log
# export DATA_FOLDER=0131
# . ./runCreateFolders.sh
# . ./runMainPipeline.sh > $LOCAL_DATASET_PATH_LOG/runMainPipeline.log

# # ./pipelineSegAndTrace_0103.sh
# # ./pipelineSegAndTrace_0113.sh #/far04
# # ./pipelineSegAndTrace_0117.sh #/far05
# # ./pipelineSegAndTrace_0120.sh #/far05
# # ./pipelineSegAndTrace_0123.sh #/far05
# ./pipelineSegAndTrace_0128.sh #/far05
# # ./pipelineSegAndTrace_0131.sh #/far04
# ./pipelineSegAndTrace_1118.sh
# ./pipelineSegAndTrace_1122.sh
# ./pipelineSegAndTrace_1230.sh
