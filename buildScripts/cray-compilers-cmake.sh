#!/bin/sh
SSD_PATH=../ssdSource
EXTRA_ARGS=$@

rm -f CMakeCache.txt
/global/u1/x/xiaofeng/cmake-2.8.5/bin/cmake \
  -D CMAKE_BUILD_TYPE:STRING=DEBUG \
  -D CMAKE_CXX_COMPILER:FILEPATH=g++ \
  -D CMAKE_Fortran_COMPILER:FILEPATH=ftn \
  -D CMAKE_Fortran_FLAGS:STRING="" \
  -D CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
  -D ALL_EXAMPLES_ENABLED:BOOL=OFF \
$EXTRA_ARGS \
$SSD_PATH
