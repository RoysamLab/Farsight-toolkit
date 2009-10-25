#!/bin/bash

# Build script for use by specific build-*.sh scripts.
# Do not call build.sh directly.

set -e
source vars.conf

if [[ -z "$PROJECT" || -z "$SRC_DIR" || -z "$BIN_DIR" ]]
then
  echo Please set PROJECT, SRC_DIR and BIN_DIR variables.
  exit 1
fi

CMAKE_CACHE=$FTK_DIR/cache/$PROJECT.cmake
OUT_FILE=$BIN_DIR/$PROJECT.txt

mkdir -p $BIN_DIR
cd $BIN_DIR
echo ======= $PROJECT: cmake =======      2>&1 | tee -a $OUT_FILE
time cmake -C $CMAKE_CACHE $SRC_DIR       2>&1 | tee -a $OUT_FILE
echo ======= $PROJECT: make =======       2>&1 | tee -a $OUT_FILE
time make                                 2>&1 | tee -a $OUT_FILE
if [ "$1" = "-install" ]
then
  echo ======= $PROJECT: install =======  2>&1 | tee -a $OUT_FILE
  time sudo make install                  2>&1 | tee -a $OUT_FILE
fi
cd -
