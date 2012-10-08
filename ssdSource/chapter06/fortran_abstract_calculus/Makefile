#The following compiler has one or more bugs preventing compilation of this code
#(Bug reports have been submitted to the compiler team) :
compile=gfortran
compile=ifort
#The following compilers correctly compile this code:
compile=xlf2003 
compile=nagfor -f2003 

executable=chaos

$(executable): main.o integrand.o lorenz.o 
	$(compile) main.o integrand.o lorenz.o -o $(executable) 

main.o: main.F90 lorenz.o 
	$(compile) -c main.F90 

lorenz.o: lorenz.F90 integrand.o
	$(compile) -c lorenz.F90

integrand.o: integrand.F90 
	$(compile) -c integrand.F90 

.PHONY : clean
clean :
	-rm -f  *.o *.mod $(executable)
