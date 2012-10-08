#Specify which Fortran code examples the chosen compiler is capable of building:
if (${ALL_EXAMPLES_ENABLED})
  add_executable(astronaut main.F90 astronaut.F90 abstract_speaker.F90 )
elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "VisualAge|XL|NAG|Intel|Cray")
  add_executable(astronaut main.F90 astronaut.F90 abstract_speaker.F90 )
elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
 message("\n***One or more reported Gnu Fortran compiler bugs prevent compilation of the Fortran ***")
 message(  "***example from Chapter 4.  Skipping build of the code in the chapter04 subdirectory.***\n")
else()
  #Build everything by default.
  add_executable(astronaut main.F90 astronaut.F90 abstract_speaker.F90 )
endif()
