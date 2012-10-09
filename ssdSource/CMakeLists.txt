#Name project and specify source languages
project(ssdSource Fortran CXX)

# Some compilers will work with earlier versions of CMake. (The Cray Compiler Environment requires CMake 2.8.5)
cmake_minimum_required(VERSION 2.8.5)

#Print an error message on an attempt to build inside the source directory tree:
IF ("${CMAKE_CURRENT_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_BINARY_DIR}")
  MESSAGE(FATAL_ERROR "ERROR! "
    "CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}"
    " == CMAKE_CURRENT_BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}"
    "\nThe ssdSource archive does not support in-source builds:\n"
    "You must now delete the CMakeCache.txt file and the CMakeFiles/ directory under"
    "the source directory for the book code archive or you will not be able to configure correctly!"
    "\nYou must now run something like:\n"
    "  $ rm -r CMakeCache.txt CMakeFiles/"
    "\n"
    "Please create a different directory and configure ssdSource under that different directory such as:\n"
    "  $ mkdir MY_BUILD\n"
    "  $ cd MY_BUILD\n"
    "  $ cmake [OPTIONS] .."
    )
ENDIF()

#Specify all the subdirectories containing codes that produces executable files:
set(chapters_to_build chapter01 chapter02 chapter03 chapter04 chapter05 chapter06 chapter07 chapter08 chapter09 chapter12)
foreach(chapter ${chapters_to_build})
  add_subdirectory(${chapter})
endforeach()

#Build a library from the directories that don't produce an executable because they contain no main program. 
#These include appendixA and chapter10.
add_subdirectory(library)

#Report untested Fortran compiler unless explicitly directed to build all examples.
if (${ALL_EXAMPLES_ENABLED})
elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "VisualAge|XL|NAG|GNU|Intel|Cray")
else()
  message("\n*** Untested Fortran compiler detected: ${CMAKE_Fortran_COMPILER_ID}.***\n")
  message("\n*** Attempting to build all examples. Please report any failures to rouson@sandia.gov.***\n")
endif()

enable_testing()

FILE(WRITE ../ssdBuild/chapter01/figure1.4/figure1.4Test.sh "./chapter01/figure1.4/figure1_4 < ../ssdSource/chapter01/figure1.4/figure1.4.in")
add_test(figure1.4Test   sh ./chapter01/figure1.4/figure1.4Test.sh)
set_tests_properties(figure1.4Test PROPERTIES PASS_REGULAR_EXPRESSION "Test passed")

add_test(figure2.4Test_1 ./chapter02/figure2.4/hello_world)
set_tests_properties(figure2.4Test_1 PROPERTIES PASS_REGULAR_EXPRESSION "Hello world!")

if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
   message("*** Skipping build of Chapter 2 figure2.4Test due to Gnu Fortran compiler bug.***\n")
else()
  add_test(figure2.4Test_2 ./chapter02/figure2.4/oo_hello_world)
  set_tests_properties(figure2.4Test_2 PROPERTIES PASS_REGULAR_EXPRESSION "Hello, world!")
endif()

FILE(WRITE ../ssdBuild/chapter02/figure2.8-11/figure2.8-11Test.sh "./chapter02/figure2.8-11/testConductor < ../ssdSource/chapter02/figure2.8-11/testConductor.in")
add_test(figure2.8-11Test sh ./chapter02/figure2.8-11/figure2.8-11Test.sh)
set_tests_properties(figure2.8-11Test PROPERTIES PASS_REGULAR_EXPRESSION "test passed")

FILE(WRITE ../ssdBuild/chapter02/figure2.14-15/figure2.14-15Test.sh "./chapter02/figure2.14-15/testFunConductor < ../ssdSource/chapter02/figure2.14-15/testFunConductor.in")
add_test(figure2.14-15Test sh ./chapter02/figure2.14-15/figure2.14-15Test.sh)
set_tests_properties(figure2.14-15Test PROPERTIES PASS_REGULAR_EXPRESSION "test passed")

if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "PGI")
   message("*** A reported PGI Fortran compiler bug prevents compilation of ***")
   message("*** chapter03/testIntegrableConductor.                     ***\n")
else()
  FILE(WRITE ../ssdBuild/chapter03/chapter03Test.sh "./chapter03/testIntegrableConductor < ../ssdSource/chapter03/testIntegrableConductor.in")
  add_test(chapter03Test sh ./chapter03/chapter03Test.sh)
  set_tests_properties(chapter03Test PROPERTIES PASS_REGULAR_EXPRESSION "test passed")
endif()

if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
   message("*** Skipping build of Chapter 4 astronautTest due to Gnu Fortran compiler bug.***\n")
else()
  add_test(astronautTest ./chapter04/astronaut)
  set_tests_properties(astronautTest PROPERTIES PASS_REGULAR_EXPRESSION "Hello, world!")
endif()

if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
   message("*** Skipping build of Chapter 5 superfluidTest due Gnu Fortran compiler bug.  ***\n")
else()
   add_test(superfluidTest ./chapter05/hermetic/superfluid)
   set_tests_properties(superfluidTest PROPERTIES PASS_REGULAR_EXPRESSION "Test passed")
endif()

add_test(chaosTest ./chapter06/fortran_abstract_calculus/chaosTest)
set_tests_properties(chaosTest PROPERTIES PASS_REGULAR_EXPRESSION "Test passed") 

if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "PGI")
   message("*** A reported Gnu Fortran compiler bug prevents compilation of ***")
   message("*** chapter07/strategy_surrogate_f2003/timed_lorenz                ***\n")
else()
   add_test(timed_lorenzTest ./chapter07/strategy_surrogate_f2003/timed_lorenz)
   set_tests_properties(timed_lorenzTest PROPERTIES PASS_REGULAR_EXPRESSION "Test passed")
endif()

add_test(weatherTest ./chapter08/puppeteer_f2003/weather)
set_tests_properties(weatherTest PROPERTIES PASS_REGULAR_EXPRESSION "Test passed")

add_test(burgersTest ./chapter09/burgers_factory_f2003/burgers)
set_tests_properties(burgersTest PROPERTIES PASS_REGULAR_EXPRESSION "Test passed")

if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  add_test(burgersCafTest ./chapter12/burgers_caf_v4/burgers_caf)
  set_tests_properties(burgersCafTest PROPERTIES PASS_REGULAR_EXPRESSION "Test passed")
endif()
