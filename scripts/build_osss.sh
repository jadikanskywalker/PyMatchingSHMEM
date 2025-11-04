#!/bin/bash

# module unload intel-oneapi-compilers
# module load gcc/14.2.0

cd ~/PyMatchingSHMEM/

CC=/mnt/DISCL/home/jadhicks/sw/el9-x86_64/osss-ucx_1.5/bin/oshcc
CXX=/mnt/DISCL/home/jadhicks/sw/el9-x86_64/osss-ucx_1.5/bin/oshcxx

export CMAKE_EXPORT_COMPILE_COMMANDS=1
cmake . -B build_osss \
  -DUSE_SHMEM=ON \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_CXX_FLAGS_DEBUG="-g -O0 -fsanitize=address" \
  -DCMAKE_C_COMPILER=$CC \
  -DCMAKE_CXX_COMPILER=$CXX \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=1

# cd build_osss
# make pymatching -j $CORES
# cd ..
