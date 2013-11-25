#Specify which Fortran code examples the chosen compiler is capable of building:
if (${ALL_EXAMPLES_ENABLED})
  add_executable(astronaut main.F90 astronaut.F90 abstract_speaker.F90 )
elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "VisualAge|XL|NAG|Intel|Cray")
  add_executable(astronaut main.F90 astronaut.F90 abstract_speaker.F90 )
elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
 message("\n***As of GCC version 4.8.1, gfortran does not yet support deferred-length   ***")
 message(  "***character components. Skipping build of the chapter04 subdirectory code. ***\n")
else()
  #Build everything by default.
  add_executable(astronaut main.F90 astronaut.F90 abstract_speaker.F90 )
endif()
