#!/bin/bash

# Script for cleaning prior build results.

set -e
source vars.conf

echo Cleaning $BUILD_DIR...
rm -rf $BUILD_DIR
