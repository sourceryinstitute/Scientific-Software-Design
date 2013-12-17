add_subdirectory(burgers_factory_cpp)

#Specify which Fortran code examples the chosen compiler is capable of building:
if (${ALL_EXAMPLES_ENABLED})
  add_subdirectory(burgers_factory_f2003)
elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "VisualAge|XL|NAG|GNU|Inel")
  add_subdirectory(burgers_factory_f2003)
elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
  message("\n***One or more reported Intel Fortran compiler bugs prevent compilation of the Fortran example from***")
  message(  "***Chapter 9.  Skipping build of the code in the chapter09/burgers_factory_f2003 subdirectory.      ***\n")
elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Cray")
  message("\n***One or more reported Cray Fortran compiler bugs prevent compilation of the Fortran example from***")
  message(  "***Chapter 9.  Skipping build of the code in the chapter09/burgers_factory_f2003 subdirectory.     ***\n")
else()
  #Build all examples by default.
  add_subdirectory(burgers_factory_f2003)
endif()
