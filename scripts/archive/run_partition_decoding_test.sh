#!/usr/bin/env bash
set -euo pipefail

# Run only the partition decoding unit test that checks virtual-boundary equivalence.
# Assumes you're in the CMake build directory (this repo root, since it uses in-tree build files).
# Usage:
#   scripts/run_partition_decoding_test.sh [--rebuild] [--build-dir DIR] [--verbose]
#
# Options:
#   --rebuild       Rebuild the pymatching_tests target before running
#   --build-dir DIR Directory containing the CMake build (default: current directory)
#   --verbose       Pass -V to ctest for verbose output

BUILD_DIR="$(pwd)"
VERBOSE=""
REBUILD="0"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --rebuild)
      REBUILD="1"; shift ;;
    --build-dir)
      BUILD_DIR="$2"; shift 2 ;;
    --verbose)
      VERBOSE="-V"; shift ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 2 ;;
  esac
done

if [[ ! -f "$BUILD_DIR/CTestTestfile.cmake" ]]; then
  echo "No CTestTestfile.cmake in $BUILD_DIR. Did you run cmake configure here?" >&2
  exit 3
fi

if [[ "$REBUILD" == "1" ]]; then
  echo "Rebuilding pymatching_tests..."
  cmake --build "$BUILD_DIR" --target pymatching_tests -j
fi

# Run just the PartitionDecoding test case
ctest --test-dir "$BUILD_DIR" -R PartitionDecoding $VERBOSE
