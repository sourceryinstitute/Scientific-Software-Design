#The following compilers have one or more bugs preventing compilation of this code
#(Bug reports have been submited to the compiler team):
compile=gfortran
compile=ifort
#The following compilers correctly compile this code:
compile=xlf2003
compile=nagfor -f2003

executable=superfluid
objects=vortex.o hermetic.o

$(executable): main.o $(objects)
	$(compile) main.o $(objects) -o $(executable) 

main.o: main.F90 vortex.o Makefile
	$(compile) -c main.F90 

vortex.o: vortex.F90 hermetic.o Makefile
	$(compile) -c vortex.F90 

hermetic.o: hermetic.F90 Makefile
	$(compile) -c hermetic.F90 

.PHONY : clean
clean :
	-rm -f  *.o *.mod $(executable) a.out core
