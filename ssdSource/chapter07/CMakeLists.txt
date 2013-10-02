add_subdirectory(strategy_cpp)

#Specify which Fortran code examples the chosen compiler is capable of building:
if (${ALL_EXAMPLES_ENABLED})
  add_subdirectory(strategy_surrogate_f2003)
elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "VisualAge|XL|NAG|GNU|Intel|PGI")
  add_subdirectory(strategy_surrogate_f2003)
else()
  #Build everything by default.
  add_subdirectory(strategy_surrogate_f2003)
endif()
