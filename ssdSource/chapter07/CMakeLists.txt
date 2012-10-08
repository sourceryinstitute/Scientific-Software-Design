add_subdirectory(strategy_cpp)

#Specify which Fortran code examples the chosen compiler is capable of building:
if (${ALL_EXAMPLES_ENABLED})
  add_subdirectory(strategy_surrogate_f2003)
elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "VisualAge|XL|NAG|GNU")
  add_subdirectory(strategy_surrogate_f2003)
elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
  message("\n***One or more reported Intel Fortran compiler bugs prevent compilation of the Fortran example from***")
  message(  "***Chapter 7.  Skipping build of the code in the chapter07/strategy_surrogate_f2003 subdirectory.   ***\n")
elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "PGI")
  message("\n***One or more reported PGI Fortran compiler bugs prevent compilation of the Fortran example from***")
  message(  "***Chapter 7.  Skipping build of the code in the chapter07/strategy_surrogate_f2003 subdirectory.  ***\n")
else()
  #Build everything by default.
  add_subdirectory(strategy_surrogate_f2003)
endif()
