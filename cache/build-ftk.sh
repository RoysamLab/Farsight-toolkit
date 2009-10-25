#!/bin/bash

# Build script for building FTK.

set -e
source vars.conf
export PROJECT=FTK
export SRC_DIR=$FTK_DIR
export BIN_DIR=$BUILD_DIR/ftk
bash build.sh
