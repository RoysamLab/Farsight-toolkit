#!/bin/bash

# Build script for building and installing ITK.

set -e
source vars.conf
export PROJECT=ITK
export SRC_DIR=$ITK_DIR
export BIN_DIR=$BUILD_DIR/itk
bash build.sh
