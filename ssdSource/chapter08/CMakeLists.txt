add_subdirectory(puppeteer_cpp)

#Specify which Fortran code examples the chosen compiler is capable of building:
if (${ALL_EXAMPLES_ENABLED})
  add_subdirectory(puppeteer_f2003)
elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "VisualAge|XL|NAG|GNU")
  add_subdirectory(puppeteer_f2003)
elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
  message("\n***One or more reported Intel Fortran compiler bugs prevent compilation of the Fortran example***")
  message(  "***from Chapter 8.  Skipping build of the code in the chapter08/puppeteer_f2003 subdirectory. ***\n")
else()
  #Build everything by default.
  add_subdirectory(puppeteer_f2003)
endif()
