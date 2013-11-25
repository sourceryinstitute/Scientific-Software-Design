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
