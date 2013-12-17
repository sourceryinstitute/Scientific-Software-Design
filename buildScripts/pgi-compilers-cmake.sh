#!/bin/sh
SSD_PATH=../ssdSource
EXTRA_ARGS=$@

rm -f CMakeCache.txt
cmake \
  -D CMAKE_BUILD_TYPE:STRING=DEBUG \
  -D CMAKE_Fortran_COMPILER:FILEPATH=pgfortran \
  -D CMAKE_Fortran_FLAGS:STRING="-Mallocatable=03" \
  -D CMAKE_VERBOSE_MAKEFILE:BOOL=ON \
  -D ALL_EXAMPLES_ENABLED:BOOL=OFF \
$EXTRA_ARGS \
$SSD_PATH
# -D CMAKE_Fortran_FLAGS:STRING="-Mallocatable=03 -lblas -L/usr/lib/" \
# -D CMAKE_CXX_COMPILER:FILEPATH=g++ \

