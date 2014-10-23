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
#GNU compiler command
compile=$(gnu)

executable=figure1.4
objects=conduction.o kind_parameters.o

$(executable): main.o $(objects) Makefile
	$(compile) main.o $(objects) -o $(executable) 

main.o: main.F90 conduction.o kind_parameters.o Makefile
	$(compile) -c main.F90 

conduction.o: conduction.F90 kind_parameters.o Makefile
	$(compile) -c conduction.F90 

kind_parameters.o: kind_parameters.F90  Makefile
	$(compile) -c kind_parameters.F90 

.PHONY : clean
clean :
	-rm -f  *.o *.mod $(executable)
