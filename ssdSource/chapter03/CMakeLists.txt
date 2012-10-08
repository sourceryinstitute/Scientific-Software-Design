#Specify which Fortran code examples the chosen compiler is capable of building:
if (${ALL_EXAMPLES_ENABLED})
  add_executable(testIntegrableConductor
    main.F90 integrable_conductor.F90 initializer.F90 field.F90 
    ../chapter02/figure2.8-11/differentiator.F90 ../chapter02/figure2.8-11/problem.F90 
    ../chapter01/figure1.4/conduction.F90 ../chapter01/figure1.4/kind_parameters.F90 
    ../appendixA/linear_solve.F90
  )
elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "VisualAge|XL")
  add_executable(testIntegrableConductor
    main.F90 integrable_conductor.F90 initializer.F90 field.F90 
    ../chapter02/figure2.8-11/differentiator.F90 ../chapter02/figure2.8-11/problem.F90 
    ../chapter01/figure1.4/conduction.F90 ../chapter01/figure1.4/kind_parameters.F90 
    ../appendixA/linear_solve.F90
  )
elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "NAG")
  add_executable(testIntegrableConductor
    main.F90 integrable_conductor.F90 initializer.F90 field.F90 
    ../chapter02/figure2.8-11/differentiator.F90 ../chapter02/figure2.8-11/problem.F90 
    ../chapter01/figure1.4/conduction.F90 ../chapter01/figure1.4/kind_parameters.F90 
    ../appendixA/linear_solve.F90
  )
elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  add_executable(testIntegrableConductor
    main.F90 integrable_conductor.F90 initializer.F90 field.F90 
    ../chapter02/figure2.8-11/differentiator.F90 ../chapter02/figure2.8-11/problem.F90 
    ../chapter01/figure1.4/conduction.F90 ../chapter01/figure1.4/kind_parameters.F90 
    ../appendixA/linear_solve.F90
  )
elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
  message("\n***One or more reported Intel Fortran compiler bugs prevent compilation of the Fortran***")
  message(  "***example from Chapter 3.  Skipping build of the code in the chapter03 subdirectory. ***\n")
elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Cray")
  message("\n***One or more reported Cray Fortran compiler bugs prevent compilation of the Fortran***")
  message(  "***example from Chapter 3.  Skipping build of the code in the chapter03 subdirectory.***\n")
elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "PGI")
  message("\n***One or more reported PGI Fortran compiler bugs prevent compilation of the Fortran***")
  message(  "***example from Chapter 3.  Skipping build of the code in the chapter03 subdirectory.***\n")
else()
  # Build everything by default
  add_executable(testIntegrableConductor
    main.F90 integrable_conductor.F90 initializer.F90 field.F90 
    ../chapter02/figure2.8-11/differentiator.F90 ../chapter02/figure2.8-11/problem.F90 
    ../chapter01/figure1.4/conduction.F90 ../chapter01/figure1.4/kind_parameters.F90 
    ../appendixA/linear_solve.F90
  )
endif()
