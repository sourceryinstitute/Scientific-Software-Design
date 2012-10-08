#This code compiles with the following compilers:
compile=xlf2003 
compile=ifort
compile=gfortran
compile=nagfor 

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
