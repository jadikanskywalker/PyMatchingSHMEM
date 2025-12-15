#!/bin/bash

cores=$(($(nproc)-2))

cd ~/PyMatchingSHMEM

cd build
make pymatching -j $cores

cd ../build_threads
make pymatching -j $cores
 
cd ..
