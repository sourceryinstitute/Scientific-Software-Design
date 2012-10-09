#Specify which Fortran code examples the chosen compiler is capable of building:
if (${ALL_EXAMPLES_ENABLED})
  add_executable(superfluid main.F90 vortex.F90 hermetic.F90 )
elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "VisualAge|XL|NAG|Cray|Intel")
  add_executable(superfluid main.F90 vortex.F90 hermetic.F90 )
elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  message("\n*** As of GCC version 4.8 20121007, gfortran does not yet support final subroutines.***")
  message(  "*** Skipping build of the code in the chapter05/hermetic subdirectory.              ***\n")
else()
  #Build everything by default.
  add_executable(superfluid main.F90 vortex.F90 hermetic.F90 )
endif()
