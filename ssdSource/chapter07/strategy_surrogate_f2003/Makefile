# The following compilers have bugs that prevent compilation of this code  
# (Bug report submitted to each vendor):
compile=gfortran      
compile=ifort         
# The following compiler compiles the code correctly:
compile=xlf2003       
compile=nagfor

objects=integrand.o lorenz.o timed_lorenz.o surrogate.o strategy.o explicit_euler.o runge_kutta_2nd.o
executable=timed_lorenz

$(executable): main.o $(objects)
	$(compile) main.o $(objects) -o $(executable) 

main.o: main.F90 $(objects)
	$(compile) -c main.F90 

surrogate.o: surrogate.F90
	$(compile) -c surrogate.F90

strategy.o: strategy.F90 surrogate.o 
	$(compile) -c strategy.F90

integrand.o: integrand.F90 surrogate.o strategy.o
	$(compile) -c integrand.F90 

lorenz.o: lorenz.F90 integrand.o
	$(compile) -c lorenz.F90 

timed_lorenz.o: timed_lorenz.F90 lorenz.o
	$(compile) -c timed_lorenz.F90 

explicit_euler.o: explicit_euler.F90 integrand.o strategy.o
	$(compile) -c explicit_euler.F90

runge_kutta_2nd.o: runge_kutta_2nd.F90 integrand.o strategy.o
	$(compile) -c runge_kutta_2nd.F90

.PHONY : clean
clean :
	-rm -f  *.o *.mod $(executable)
