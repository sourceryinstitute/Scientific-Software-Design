#!/bin/sh
SSD_PATH=../ssdSource
EXTRA_ARGS=$@

rm -f CMakeCache.txt
cmake \
  -D CMAKE_BUILD_TYPE:STRING=DEBUG \
  -D CMAKE_CXX_COMPILER:FILEPATH=g++ \
  -D CMAKE_Fortran_COMPILER:FILEPATH=gfortran \
  -D CMAKE_Fortran_FLAGS:STRING="-lblas -llapack -fcoarray=single" \
  -D CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
  -D ALL_EXAMPLES_ENABLED:BOOL=OFF \
$EXTRA_ARGS \
$SSD_PATH
