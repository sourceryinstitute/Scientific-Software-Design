# Choose desired compiler (options: gnu, intel, cray, tau_intel):
compiler=gnu
executable=burgers_caf
mofo_root=/Users/rouson/mofo
fireworks=$(mofo_root)/Examples/serial_examples/src/chapter01/fireworks
library_root=./library/
include_root=./include-files/

ifeq ($(compiler),gnu)
  opt=-L /opt/opencoarrays/lib
  compile=mpif90 -fcoarray=lib
  lib=-lcaf_mpi
else
  ifeq ($(compiler),intel)
    compile=ifort -coarray=shared -standard-semantics -O3 -coarray-num-images=2
  else
    ifeq ($(compiler),cray)
      compile=ftn -ew -h caf
    else
      ifeq ($(compiler),tau_intel)
        compile=tauf90 -tau:serial,icpc,pdt  -coarray=shared -standard-semantics -O3 -coarray-num-images=4 -DTAU=1
      else
        ifeq ($(compiler),tau_aciss)
          tau_aciss=tau_f90.sh -v -Bdynamic -coarray=shared -standard-semantics -O3
        else
          ifeq ($(compiler),tau_cray)
            tau_hopper=tau_f90.sh -optCompInst -ew -h caf
          endif
        endif
      endif
    endif
  endif
endif

objects = local_field.o global_field.o kind_parameters.o co_object_interface.o object_interface.o ForTrilinos_assertion_utility.o ForTrilinos_error.o math_constants.o burgers.o

$(executable): main.o $(objects) Makefile
	$(compile) $(opt)  main.o $(objects) -o $(executable)  $(lib)


main.o: main.F90 $(objects) Makefile
	$(compile) -c main.F90

global_field.o: global_field.F90 local_field.o kind_parameters.o co_object_interface.o ForTrilinos_assertion_utility.o ForTrilinos_error.o preprocessor_definitions Makefile
	$(compile) -c global_field.F90

local_field.o: local_field.F90 kind_parameters.o object_interface.o ForTrilinos_assertion_utility.o ForTrilinos_error.o Makefile
	$(compile) -c local_field.F90

burgers.o: burgers.F90 kind_parameters.o object_interface.o ForTrilinos_assertion_utility.o ForTrilinos_error.o math_constants.o Makefile
	$(compile) -c burgers.F90

math_constants.o: math_constants.o kind_parameters.o Makefile
	$(compile) -c math_constants.F90

kind_parameters.o: kind_parameters.F90 Makefile
	$(compile) -c kind_parameters.F90

object_interface.o: $(library_root)/object_interface.F90 preprocessor_definitions Makefile 
	$(compile) -c $(library_root)/object_interface.F90

co_object_interface.o: $(library_root)/co_object_interface.F90 Makefile
	$(compile) -c $(library_root)/co_object_interface.F90

ForTrilinos_assertion_utility.o: $(library_root)/ForTrilinos_assertion_utility.F90 Makefile
	$(compile) -c $(library_root)/ForTrilinos_assertion_utility.F90

ForTrilinos_error.o: $(library_root)/ForTrilinos_error.F90 Makefile
	$(compile) -c $(library_root)/ForTrilinos_error.F90

.PHONY : clean
clean :
	-rm -f  *.o *.mod $(executable) core fort.* compiler_capabilities.txt $(library_root)/compiler_capabilities.txt profile* MULT* *.trc *.edf

preprocessor_definitions: 
ifeq ($(compiler),gnu)
	cp $(include_root)/gfortran_capabilities.txt compiler_capabilities.txt
	cp $(include_root)/gfortran_capabilities.txt $(library_root)/compiler_capabilities.txt
else
  ifeq ($(compiler),intel)
	cp $(include_root)/intel_capabilities.txt compiler_capabilities.txt
	cp $(include_root)/intel_capabilities.txt $(library_root)/compiler_capabilities.txt
  else
    ifeq ($(compiler),cray)
	cp $(include_root)/cray_capabilities.txt compiler_capabilities.txt
	cp $(include_root)/cray_capabilities.txt $(library_root)/compiler_capabilities.txt
  else
    ifeq ($(compiler),tau_cray)
	cp $(include_root)/tau_cray_capabilities.txt compiler_capabilities.txt
	cp $(include_root)/tau_cray_capabilities.txt $(library_root)/compiler_capabilities.txt
  else
    ifeq ($(compiler),tau_intel)
	cp $(include_root)/tau_intel_capabilities.txt compiler_capabilities.txt
	cp $(include_root)/tau_intel_capabilities.txt $(library_root)/compiler_capabilities.txt
    endif
  endif
  endif
  endif
endif
