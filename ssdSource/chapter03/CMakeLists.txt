#Specify which Fortran code examples the chosen compiler is capable of building:
if (${ALL_EXAMPLES_ENABLED})
  add_executable(testIntegrableConductor
    main.F90 integrable_conductor.F90 initializer.F90 field.F90 
    ../chapter02/figure2.8-11/differentiator.F90 ../chapter02/figure2.8-11/problem.F90 
    ../chapter01/figure1.4/conduction.F90 ../chapter01/figure1.4/kind_parameters.F90 
    ../appendixA/linear_solve.F90
  )
elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "VisualAge|XL|NAG|GNU|PGI|Cray")
  add_executable(testIntegrableConductor
    main.F90 integrable_conductor.F90 initializer.F90 field.F90 
    ../chapter02/figure2.8-11/differentiator.F90 ../chapter02/figure2.8-11/problem.F90 
    ../chapter01/figure1.4/conduction.F90 ../chapter01/figure1.4/kind_parameters.F90 
    ../appendixA/linear_solve.F90
  )
elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
  message("\n*** A reported Intel compiler bug prevents compilation of the chapter03 ***")
  message(  "*** code:  Intel Premier Support -- Issue number 6000036099.            ***\n")
endif()
