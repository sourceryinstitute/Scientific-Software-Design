The Scientific Software Design Repository
=========================================

Overview
--------
This repository contains modern Fortran  and C++ codes first published in the textbook 
Scientific Software Design: The Object-Oriented Way by Damian Rouson, Jim Xia, and
Xiaofeng Xu (Cambridge University Press, 2011).  Updated versions of this repository
can be obtained via git as follows:

  git clone https://your-user-name@bitbucket.org/sourceryinstitute/scientific-software-design.git

A less-frequently updated version of this code archive also resides at 
http://wwww.cambridge.org/Rouson via the "Resources" link.

Organization
------------
The codes are in the ssdSource directory and organized into subdirectories named 
after the chapter in which the code can be found.  The ssdSource directoy holds a 
CMake project (http://www.cmake.org). The buildScripts directory contains CMake 
commands for each supported compiler.  The ssdBuild directory, if present, can be 
used for building all of the executables in an automated fashion.  (If ssdBuild 
does not exist, create it or a directory with the name of your choice at the same 
level as the buildScripts directory before proceeding to the "Building" section 
below.)

Building
--------
If you are building inside the HPCLinux distribution available at 
http://www.hpclinux.org, please enter the command "source ./intel.sh" inside the
scientific-software-design directory before continuing.

The preferred method for compiling the code in this repository relies upon the 
automated, open-source CMake build system.  CMake is available for multiple operating
systems at http://www.cmake.org.  CMake can also be installed automatically via 
package management software such as Macports for OS X or yum or apt-get for Linux.  
After installing CMake, set your current directory to this archive's ssdBuild 
subdirectory and run the build script corresponding to your compiler (after editing 
the script to update its details to match your configuration) using a command such as

cd ssdBuild
../buildScript/gnu-compilers-cmake.sh
make
make test

The above commands would use the GNU Fortran, C and C++ compilers to build the 
executable files and would test each executable and report on the success of each 
test.  Alternatively, one can build just a given chapter's code.  For example, to
build the Chapter 3 examples only, enter the following commands:

../buildScript/gnu-compilers-cmake.sh
cd chapter03
make

The aforementioned editing of the script must include providing a path to the LAPACK
library (http://www.netlib.org/lapack/), which is required to build code in Chapter 9.

Most folders also contain a Makefile for building with Unix-like make utilities.
Because it requires a great deal more effort to maintain Makefiles relative to the 
effort of maintaining CMake files, we make no promise that the Makefiles will always
be maintained and up-to-date -- hence the recommendation that CMake is the preferred
method for building the archive.

Most subdirectories in this archive also contain Makefiles that can be used to compile
the code on platforms with the "make" utility 
(http://en.wikipedia.org/wiki/Make_(software)).  We deprecate the use of these 
Makefiles because of the amount of effort that would be required to keep them 
up-to-date and the number of files users must customize to build the full archive.  

Compiler Prerequisites
----------------------
This archive makes extensive use of the explicit support for object-oriented 
programming in the Fortran 2003 standard (and Fortran 2008 in one example in Chapter
12).  As of November 2013, six compilers nominally support all of the Fortran 2003
features employed in this code archive:
  1. IBM XL (xlf),
  2. The Cray Compiler Environment (cce),
  3. Portland Group (pgfortran),
  4. Numerical Algorithms Group (nagfor)
  5. Intel (ifort)
  6. GNU (gfortran)
Most compilers not listed above do not support a sufficient number of Fortran 2003
features to compile most code in the book.  If your compiler is listed above but you
are unable to compile the code in ssdSource, please update to the latest version of
your compiler and retry before reporting the bug to damian@sourceryinstitute.org.

Additionally, the Cray and Intel compilers and the GNU 5.0 compiler fully support the 
Fortran 2008 parallel programming features needed to build the coarray example in 
Chapter 12.  Producing parallel executables with GNU Fortran requires building and
linking to an OpenCoarrays (http://opencoarrays.org) transport layer.
