# The following compilers do not compile this code correctly
# (Bug reports have been submitted to the compiler teams):
compile=gfortran
compile=ifort
# The following compilers compile this code correctly.
compile=xlf2003 
compile=nagfor -f2003 

executable=weather

$(executable): main.o atmosphere.o air.o cloud.o ground.o global_parameters.o
	$(compile) main.o atmosphere.o air.o cloud.o ground.o -o $(executable) 

main.o: main.F90 atmosphere.o global_parameters.o
	$(compile) -c main.F90 

atmosphere.o: atmosphere.F90 air.o ground.o cloud.o global_parameters.o
	$(compile) -c atmosphere.F90 

air.o: air.F90 global_parameters.o
	$(compile) -c air.F90 

cloud.o: cloud.F90 global_parameters.o
	$(compile) -c cloud.F90 

ground.o: ground.F90 global_parameters.o
	$(compile) -c ground.F90 

global_parameters.o: global_parameters.F90 
	$(compile) -c global_parameters.F90 

.PHONY : clean
clean :
	-rm -f  *.o *.mod $(executable)
