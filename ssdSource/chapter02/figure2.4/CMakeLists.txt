add_executable(hello_world hello_world.F90)

if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
   message("*** As of GCC version 4.8.1, gfortran does not yet support using deferred-length     ***")
   message("*** character components.  Skipping build of chapter02/figure2.4/oo_hello_world.F90. ***\n")
else()
   add_executable(oo_hello_world oo_hello_world.F90)
endif()
