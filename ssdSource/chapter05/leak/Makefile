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

executable=leak
objects=leak.o

$(executable): main.o $(objects)
	$(compile) main.o $(objects) -o $(executable) 

main.o: main.F90 leak.o Makefile
	$(compile) -c main.F90 

leak.o: leak.F90 Makefile
	$(compile) -c leak.F90 

.PHONY : clean
clean :
	-rm -f  *.o *.mod $(executable) a.out core
