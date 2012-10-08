This code archive derives from the examples in the textbook Scientific Software Design: The Object-Oriented Way by Damian Rouson, Jim Xia, and Xiaofeng Xu (Cambridge University Press, 2011).  The preferred method for compiling the code in this directory relies upon the automated, open-source CMake build system (http://www.cmake.org).  CMake can also be installed automatically via package management software such as Macports for OS X or yum for Linux.  After installing CMake, set your current directory to this archive's ssdBuild subdirectory and run the build script corresponding to your compiler (after editing the script to update its details to match your configuration) using a command such as

../buildScript/gnu-compilers-cmake.sh
make

The above command would use the GNU Fortran and C++ compilers to build the executable files and the library file. The aforementioned editing of the script must include providing a path to the LAPACK library (http://www.netlib.org/lapack/), which is required to build code in Chapter 9.  Alternatively, one can build just a given chapter's code.  For example, to build the 

../buildScript/gnu-compilers-cmake.sh
cd chapter06/fortran_abstract_calculus 
make

In this case, one can build every supplying a LAPACK build in the build script except for the chapter09 subdirectory.

Most subdirectories in this archive also contain Makefiles that can be used to compile the code on platforms with the "make" utility (http://en.wikipedia.org/wiki/Make_(software)).  Each Makefile indicates which compilers correctly compiled the code as of August 2011.  We deprecate the use of these Makefiles because of the amount of effort that would be required to keep them up-to-date and the number of files users must customize to build the full archive.  

This archive makes extensive use of the explicit support for object-oriented programming in the Fortran 2003 standard (and Fortran 2008 in one example in Chapter 12).  As of August 2011, three compilers nominally support the complete Fortran 2003 standard:
  1. IBM XL (xlf)
  2. The Cray Compiler Environment (cce)
  3. Portland Group (pgfortran)
while three compilers nominally support all of the Fortran 2003 features used in this code archive:
  4. Numerical Algorithms Group (nagfor)
  5. Intel (ifort)
and one compiler supports all Fortran 2003 features used in this archive except two:
  6. GNU (gfortran)
Most compilers not listed above do not support a sufficient number of Fortran 2003 features to compile most code in the book.

For every bug or missing feature of which we are aware, we have submitted a problem report to the appropriate compiler team and adjusted the CMake scripts to disable the building of the corresponding code.  We would greatly appreciate it if users of this code would report any new problems directly to the appropriate compiler team as wel as notifying the first author at damian@rouson.net.

We encourage readers interested in compiling with an open-source compiler to contribute code to the gfortran project and let the gfortran developers know about your interest in resolving the following bug reports by e-mailing fortran@gcc.gnu.org:
* Bug 47545: This prevents compiling with deferred-length components. A workaround is chapter02/figure2.4/gfortran_oo_hello_world.F90.
* Bug 45170 (comment 9): This presumably precludes returning character variables with allocatable lengths.
* Bug 37336: This prevents use of final subroutines. Chapter 5 describes a workaround of sorts.
* Bug 18918: This prevents use of Fortran 2008 coarrays in chapter12/burgers_caf_v4. No workaround is provided. (This bug is marked as fixed but gfortran cannot yet compile the burgers_caf_v4 code.)

Any updates to this archive will be accessible via the "Resources Available" link at http://wwww.cambridge.org/Rouson.
