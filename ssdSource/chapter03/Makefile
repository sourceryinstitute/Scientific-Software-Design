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


executable=testIntegrableConductor
objects=conduction.o differentiator.o problem.o kind_parameters.o integrable_conductor.o initializer.o field.o linear_solve.o

$(executable): main.o $(objects)
	$(compile) main.o $(objects) -o $(executable) 

main.o: main.F90 kind_parameters.o integrable_conductor.o initializer.o
	$(compile) -c main.F90 

initializer.o: initializer.F90 kind_parameters.o 
	$(compile) -c initializer.F90 

integrable_conductor.o: integrable_conductor.F90 field.o problem.o kind_parameters.o 
	$(compile) -c integrable_conductor.F90 

field.o: field.F90 kind_parameters.o differentiator.o problem.o linear_solve.o
	$(compile) -c field.F90 

linear_solve.o: ../appendixA/linear_solve.F90 kind_parameters.o 
	$(compile) -c ../appendixA/linear_solve.F90 

differentiator.o: ../chapter02/figure2.8-11/differentiator.F90 conduction.o problem.o kind_parameters.o
	$(compile) -c ../chapter02/figure2.8-11/differentiator.F90 

problem.o: ../chapter02/figure2.8-11/problem.F90 conduction.o kind_parameters.o
	$(compile) -c ../chapter02/figure2.8-11/problem.F90 

conduction.o: ../chapter01/figure1.4/conduction.F90 kind_parameters.o
	$(compile) -c ../chapter01/figure1.4/conduction.F90 

kind_parameters.o: ../chapter01/figure1.4/kind_parameters.F90 
	$(compile) -c ../chapter01/figure1.4/kind_parameters.F90 

.PHONY : clean
clean :
	-rm -f  *.o *.mod $(executable)
	-rm -f  ../chapter02/figure2.8-11/*.o ../chapter02/figure2.8-11/*.mod 
	-rm -f  ../chapter01/figure1.4/*.o ../chapter01/figure1.4/*.mod 
