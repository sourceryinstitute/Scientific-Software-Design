if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  message("\n")
  message("*** As of GCC 4.8 20121007, gfortran does not support the deferred-length character ***")
  message("*** components required by chapter10/assertions/assertion_utility.F90               ***\n")

  message("*** Nor does gfortran yet support the final subroutines required by                       ***")
  message("*** chapter11/reference_count/universal.F90 and chapter11/reference_count/ref_counter.F90 ***")
  message("*** These codes will not be added to the libssd.a library.                                ***\n")
  add_library(ssd 
     ../appendixA/linear_solve.F90 
     ../chapter01/figure1.4/kind_parameters.F90
     ../chapter11/reference_count/hermetic.F90
  ) 
else()
  add_library(ssd 
     ../chapter10/assertions/assertion_utility.F90 
     ../appendixA/linear_solve.F90 
     ../chapter01/figure1.4/kind_parameters.F90
     ../chapter11/reference_count/hermetic.F90
     ../chapter11/reference_count/universal.F90
     ../chapter11/reference_count/ref_counter.F90
  ) 
endif()
