#!/bin/bash

# Build script for building all components (VXL, ITK, VTK and FTK).

set -e
bash build-vxl.sh
bash build-itk.sh
bash build-vtk.sh
bash build-ftk.sh
