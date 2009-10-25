#!/bin/bash

# Build script for building and installing VXL.

set -e
source vars.conf
export PROJECT=VXL
export SRC_DIR=$VXL_DIR
export BIN_DIR=$BUILD_DIR/vxl
bash build.sh
