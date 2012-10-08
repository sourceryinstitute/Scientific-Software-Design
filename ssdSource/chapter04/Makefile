#This following compiler has one or more bugs that prevent compilation of this code
#(Bug reports have been submitted to the compiler team):
compile=gfortran
#This code works with the following compilers:
compile=xlf2003
compile=ifort
compile=nagfor -f2003

executable=astronaut
objects=astronaut.o abstract_speaker.o

$(executable): main.o $(objects)
	$(compile) main.o $(objects) -o $(executable) 

main.o: main.F90 astronaut.o Makefile
	$(compile) -c main.F90 

astronaut.o: astronaut.F90 abstract_speaker.o Makefile
	$(compile) -c astronaut.F90 

abstract_speaker.o: abstract_speaker.F90 Makefile
	$(compile) -c abstract_speaker.F90 

.PHONY : clean
clean :
	-rm -f  *.o *.mod $(executable)
