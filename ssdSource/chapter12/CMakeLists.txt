#Specify which Fortran code examples the chosen compiler is capable of building:
if (${ALL_EXAMPLES_ENABLED})
  add_subdirectory(burgers_caf_v5)
elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Cray|GNU|Intel")
  add_subdirectory(burgers_caf_v5)
else ()
  message(  "***Chapter 12.  Skipping build of the code in the chapter12/burgers_caf_v5 subdirectory.      ***\n")
endif()
