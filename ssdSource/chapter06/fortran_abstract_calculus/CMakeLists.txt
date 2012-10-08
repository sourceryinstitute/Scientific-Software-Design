if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  message("\n***One or more reported GNU Fortran compiler bugs prevent compilation of the Fortran example from***")
  message(  "***Chapter 6.  Skipping build of the code in the chapter06/fortran_abstract_calculus subdirectory.  ***\n")
else()
  #Build everything by default.
  add_executable(chaos main.F90  lorenz.F90 integrand.F90)
endif()
