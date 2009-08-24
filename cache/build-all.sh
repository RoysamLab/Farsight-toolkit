#!/bin/bash

# NB: Change these paths to point to your copy
# of the source code for the relevant packages.
FTK_DIR=~/svn/farsight
ITK_DIR=~/cvs/Insight
VTK_DIR=~/src/VTK
VXL_DIR=~/src/vxl-1.12.0
BOOST_DIR=/usr/local/include/boost-1_38

bash build-vxl.sh
bash build-itk.sh
bash build-vtk.sh
bash build-ftk.sh
