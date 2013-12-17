#!/bin/sh
SSD_PATH=../ssdSource
EXTRA_ARGS=$@

rm -f CMakeCache.txt
cmake \
  -D CMAKE_Fortran_COMPILER:FILEPATH=ifort \
  -D CMAKE_Fortran_FLAGS:STRING="-standard-semantics -coarray=shared -coarray-num-images=1" \
  -D CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
  -D ALL_EXAMPLES_ENABLED:BOOL=OFF \
$EXTRA_ARGS \
$SSD_PATH
#  -D CMAKE_CXX_COMPILER:FILEPATH=iCC \
#  -D CMAKE_BUILD_TYPE:STRING=DEBUG \
#  -D CMAKE_Fortran_FLAGS:STRING="-standard-semantics -coarray=shared -coarray-num-images=1 -lblas -llapack -L/usr/lib/libblas:/usr/lib/liblapack" \
