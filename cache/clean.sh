#!/bin/bash

# Script for cleaning prior build results.

set -e
source vars.conf

rm -rf $BUILD_DIR
