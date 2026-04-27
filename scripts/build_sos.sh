#!/bin/bash

# Default to debug build if no argument is provided
BUILD_TYPE=${1:-debug}
FLAG_TYPE=${1:-fast}

cd ~/PyMatchingSHMEM/

# rm -rf build_sos

CC=/mnt/DISCL/home/jadhicks/sw/el9-x86_64/sos_1.5_scalable/bin/oshcc
CXX=/mnt/DISCL/home/jadhicks/sw/el9-x86_64/sos_1.5_scalable/bin/oshc++

# if [ "$FLAG_TYPE" = "fast" ]; then
#     RELEASE_FLAGS="-O2 -fno-omit-frame-pointer -g"
# else
#     RELEASE_FLAGS="-O3 -DNDEBUG"
# fi

export CMAKE_EXPORT_COMPILE_COMMANDS=1
if [ "$BUILD_TYPE" = "release" ]; then
    rm -rf build_sos
    echo "Configuring Release build..."
    cmake . -B build_sos \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_CXX_FLAGS_RELEASE="$RELEASE_FLAGS" \
        -DCMAKE_EXPORT_COMPILE_COMMANDS=1 \
        -DCMAKE_C_COMPILER="$CC" \
        -DCMAKE_CXX_COMPILER="$CXX" \
        -DUSE_THREADS=ON \
        -DUSE_SHMEM=ON
    cd build_sos|| exit
elif [ "$BUILD_TYPE" = "debug" ]; then
    rm -rf build_sos
    echo "Configuring Debug build..."
    cmake . -B build_sos \
        -DCMAKE_BUILD_TYPE=Debug \
        -DCMAKE_CXX_FLAGS_DEBUG="-g -O3" \
        -DCMAKE_EXPORT_COMPILE_COMMANDS=1 \
        -DCMAKE_C_COMPILER=$CC \
        -DCMAKE_CXX_COMPILER=$CXX \
        -DUSE_THREADS=ON \
        -DUSE_SHMEM=ON
    cd build_sos || exit
elif [ "$BUILD_TYPE" = "profile" ]; then
    rm -rf build_sos_profile
    echo "Configuring Debug+Profile build..."
    cmake . -B build_sos_profile \
        -DCMAKE_BUILD_TYPE=Debug \
        -DCMAKE_CXX_FLAGS_DEBUG="-g -O3" \
        -DCMAKE_EXPORT_COMPILE_COMMANDS=1 \
        -DCMAKE_C_COMPILER=$CC \
        -DCMAKE_CXX_COMPILER=$CXX \
        -DCMAKE_C_COMPILER_LAUNCHER="scorep;--user;--mpp=none;--thread=none" \
        -DCMAKE_CXX_COMPILER_LAUNCHER="scorep;--user;--mpp=none;--thread=none" \
        -DCMAKE_C_LINKER_LAUNCHER="scorep;--user;--mpp=none;--thread=none" \
        -DCMAKE_CXX_LINKER_LAUNCHER="scorep;--user;--mpp=none;--thread=none" \
        -DUSE_THREADS=ON \
        -DUSE_SHMEM=ON
    cd build_sos_profile || exit
else
    echo "Invalid build type: $BUILD_TYPE. Use 'debug' or 'release'."
    exit 1
fi


cores=$(($(nproc)-1))
make pymatching -j "$cores"
cd ..
