#!/bin/bash

# module unload intel-oneapi-compilers
# module load gcc/14.2.0

cd ~/PyMatchingSHMEM/

CC=$(which gcc)
CXX=$(which g++)

export CMAKE_EXPORT_COMPILE_COMMANDS=1
cmake . -B build_threads \
  -DUSE_THREADS=ON \
  -DUSE_SHMEM=OFF \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_CXX_FLAGS_DEBUG="-g -O0 -fsanitize=address" \
  -DCMAKE_C_COMPILER=$CC \
  -DCMAKE_CXX_COMPILER=$CXX \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=1

# cd build_osss
# make pymatching -j $CORES
# cd ..
