executable=burgers_caf
#Cray compiler command
cray=ftn -ew -h caf
#Intel compiler command
intel=ifort -standard-semantics -coarray=shared -coarray-num-images=4
#IBM compiler command
ibm=xlf2003
#NAG compiler command
nag=nagfor -f2008 -fpp -C=all
#GNU compiler command
gnu=gfortran -std=f2008 -fcheck=all -fcoarray=single -Wall -Wno-unused-dummy-argument
#Portland Group compiler command
pg=pgfortran -Mallocatable=03
#GNU compiler command
compile=$(intel)

objects = field.o periodic_2nd_order.o kind_parameters.o 

$(executable): main.o $(objects) Makefile
	$(compile)  main.o $(objects) -o $(executable) 

main.o: main.F90 $(objects) Makefile
	$(compile) -c main.F90

periodic_2nd_order.o: periodic_2nd_order.F90 field.o kind_parameters.o Makefile
	$(compile) -c periodic_2nd_order.F90

field.o: field.F90 kind_parameters.o Makefile
	$(compile) -c field.F90

kind_parameters.o: kind_parameters.F90 Makefile
	$(compile) -c kind_parameters.F90

.PHONY : clean
clean :
	-rm -f  *.o *.mod $(executable) core
