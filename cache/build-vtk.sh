#!/bin/bash

# Build script for building and installing VTK.

set -e
source vars.conf
export PROJECT=VTK
export SRC_DIR=$VTK_DIR
export BIN_DIR=$BUILD_DIR/vtk
bash build.sh

# HACK - Fix for broken VTK install target that fails to copy certain files.
#VTK_INCLUDE_DIR=/usr/local/include/vtk-5.4
#sudo cp $VTK_DIR/GUISupport/Qt/vtkQtBarChartView.h $VTK_INCLUDE_DIR
#sudo cp $VTK_DIR/GUISupport/Qt/vtkQtChartViewBase.h $VTK_INCLUDE_DIR
