if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  message("\n***One or more Gnu Fortran Compiler bugs prevent compilation of ***")
  message("***   chapter10/assertions/assertion_utility.F90                ***")
  message("***   chapter11/reference_count/universal.F90                   ***")
  message("***   chapter11/reference_count/ref_counter.F90                 ***")
  message("***These codes will not be added to the libssd.a library.         ***\n")
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
