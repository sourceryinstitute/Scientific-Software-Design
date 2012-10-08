#The following compilers correctly compile this code:
compile=xlf2003
compile=ifort
compile=nagfor -f2003
compile=gfortran

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
