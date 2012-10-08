#The following compiler compiles this code but generates a runtime error:
compile=ifort
#The following compilers correctly compiler this code:
compile=xlf2003
compile=nagfor -f2003
compile=gfortran

executable=testConductor
objects=conduction.o differentiator.o problem.o conductor.o kind_parameters.o

$(executable): main.o $(objects)
	$(compile) main.o $(objects) -o $(executable) 

main.o: main.F90 differentiator.o problem.o conductor.o kind_parameters.o
	$(compile) -c main.F90 

conductor.o: conductor.F90 conduction.o kind_parameters.o differentiator.o problem.o
	$(compile) -c conductor.F90 

differentiator.o: differentiator.F90 conduction.o problem.o kind_parameters.o
	$(compile) -c differentiator.F90 

problem.o: problem.F90 conduction.o kind_parameters.o
	$(compile) -c problem.F90 

conduction.o: ../../chapter01/figure1.4/conduction.F90 kind_parameters.o
	$(compile) -c ../../chapter01/figure1.4/conduction.F90 

kind_parameters.o: ../../chapter01/figure1.4/kind_parameters.F90 
	$(compile) -c ../../chapter01/figure1.4/kind_parameters.F90 

.PHONY : clean
clean :
	-rm -f  *.o *.mod $(executable)
