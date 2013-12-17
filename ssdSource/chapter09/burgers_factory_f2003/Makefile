#Cray compiler command
cray=ftn -ew -h
#Intel compiler command
intel=ifort -standard-semantics 
#IBM compiler command
ibm=xlf2003
#NAG compiler command
nag=nagfor -f2008 -fpp -C=all
#GNU compiler command
gnu=gfortran -std=f2008
#Portland Group compiler command
pg=pgfortran -Mallocatable=03
#Edit the following line to select the compiler
compile=$(gnu)

ifeq ($(compile),$(ibm))
  library = -lessl
else
  ifeq ($(compile),$(cray))
    library = 
  else
    library = -lblas -llapack -L/usr/lib
  endif
endif

executable=burgers
objects = field.o field_factory.o periodic_6th_factory.o periodic_6th_order.o matrix.o kind_parameters.o 

$(executable): main.o $(objects) Makefile
	$(compile)  main.o $(objects) -o $(executable) $(library)

main.o: main.F90 $(objects) Makefile
	$(compile) -c main.F90

periodic_6th_factory.o: periodic_6th_factory.F90 periodic_6th_order.o field_factory.o kind_parameters.o Makefile
	$(compile) -c periodic_6th_factory.F90

periodic_6th_order.o: periodic_6th_order.F90 field.o matrix.o kind_parameters.o Makefile
	$(compile) -c periodic_6th_order.F90

field_factory.o: field_factory.F90 field.o kind_parameters.o Makefile
	$(compile) -c field_factory.F90

field.o: field.F90 kind_parameters.o Makefile
	$(compile) -c field.F90

matrix.o: matrix.F90 kind_parameters.o Makefile
	$(compile) -c matrix.F90

kind_parameters.o: ../../chapter01/figure1.4/kind_parameters.F90 Makefile
	$(compile) -c ../../chapter01/figure1.4/kind_parameters.F90

.PHONY : clean
clean :
	-rm -f  *.o *.mod $(executable) core
