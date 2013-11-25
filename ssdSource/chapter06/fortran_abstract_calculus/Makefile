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
