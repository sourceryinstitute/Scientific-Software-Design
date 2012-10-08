#!/bin/sh
SSD_PATH=../ssdSource
EXTRA_ARGS=$@

rm -f CMakeCache.txt
cmake \
  -D CMAKE_BUILD_TYPE:STRING=DEBUG \
  -D CMAKE_CXX_COMPILER:FILEPATH=g++ \
  -D CMAKE_Fortran_COMPILER:FILEPATH=nagfor \
  -D CMAKE_Fortran_FLAGS:STRING="-lblas -llapack -L/Users/rouson/bin/lapack-3.4.1/" \
  -D CMAKE_VERBOSE_MAKEFILE:BOOL=ON \
  -D ALL_EXAMPLES_ENABLED:BOOL=OFF \
$EXTRA_ARGS \
$SSD_PATH
