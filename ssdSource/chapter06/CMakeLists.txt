add_subdirectory(semidiscrete_cpp)

#Specify which Fortran code examples the chosen compiler is capable of building:
if (${ALL_EXAMPLES_ENABLED})
  add_subdirectory(fortran_abstract_calculus)
elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "VisualAge|XL|NAG|GNU|PGI|Cray|Intel")
  add_subdirectory(fortran_abstract_calculus)
else()
  #Build everything by default.
  add_subdirectory(fortran_abstract_calculus)
endif()
