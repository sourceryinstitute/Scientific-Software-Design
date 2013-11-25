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

executable=testFunConductor
objects=conduction.o differentiator.o problem.o conductor.o kind_parameters.o fun_conductor.o

$(executable): main.o $(objects)
	$(compile) main.o $(objects) -o $(executable) 

main.o: main.F90 kind_parameters.o fun_conductor.o
	$(compile) -c main.F90 

fun_conductor.o: fun_conductor.F90 conductor.o differentiator.o problem.o kind_parameters.o
	$(compile) -c fun_conductor.F90 

conductor.o: ../figure2.8-11/conductor.F90 conduction.o kind_parameters.o differentiator.o
	$(compile) -c ../figure2.8-11/conductor.F90 

differentiator.o: ../figure2.8-11/differentiator.F90 conduction.o problem.o kind_parameters.o
	$(compile) -c ../figure2.8-11/differentiator.F90 

problem.o: ../figure2.8-11/problem.F90 conduction.o kind_parameters.o
	$(compile) -c ../figure2.8-11/problem.F90 

conduction.o: ../../chapter01/figure1.4/conduction.F90 kind_parameters.o
	$(compile) -c ../../chapter01/figure1.4/conduction.F90 

kind_parameters.o: ../../chapter01/figure1.4/kind_parameters.F90 
	$(compile) -c ../../chapter01/figure1.4/kind_parameters.F90 

.PHONY : clean
clean :
	-rm -f  *.o *.mod $(executable)
