#!/bin/bash

cd ~/PyMatchingSHMEM

export ASAN_OPTIONS="detect_leaks=0:halt_on_error=0"
echo $ASAN_OPTIONS
cmake . -DCMAKE_BUILD_TYPE=Debug

make clean
make pymatching -j $CORES