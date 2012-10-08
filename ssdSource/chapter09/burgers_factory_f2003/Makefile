#The following compilers have bugs that prevent compilation of this code
#(Bug reports have been submitted to the compiler teams):
compile=ifort
compile=gfortran
compile=nagfor -f2003
#The following compilers compile this code correctly:
compile=xlf2003 
compile=nagfor

executable=burgers
objects = field.o field_factory.o periodic_6th_factory.o periodic_6th_order.o matrix.o kind_parameters.o 
#Use this line if you have built LAPACK/BLAS from source obtained from http://www.netlib.org.
#library = -lblas -L/path-to-lapack/lapack-3.3.1/
library = -lblas -llapack -L/Users/rouson/bin/lapack-3.4.1/
#Use this line if you will like to the LAPACK/BLAS optimized for IBM POWER architecture platforms:
#library = -lessl

$(executable): main.o $(objects) Makefile
	$(compile)  main.o $(objects) -o $(executable) $(library)

main.o: main.F90 $(objects) Makefile
	$(compile) -c main.F90

periodic_6th_factory.o: periodic_6th_factory.F90 periodic_6th_order.o field_factory.o kind_parameters.o Makefile
	$(compile) -c periodic_6th_factory.F90

periodic_6th_order.o: periodic_6th_order.F90 field.o matrix.o kind_parameters.o Makefile
	$(compile) -c periodic_6th_order.F90

field_factory.o: field_factory.F90 field.o kind_parameters.o Makefile
	$(compile) -c field_factory.F90

field.o: field.F90 kind_parameters.o Makefile
	$(compile) -c field.F90

matrix.o: matrix.F90 kind_parameters.o Makefile
	$(compile) -c matrix.F90

kind_parameters.o: ../../chapter01/figure1.4/kind_parameters.F90 Makefile
	$(compile) -c ../../chapter01/figure1.4/kind_parameters.F90

.PHONY : clean
clean :
	-rm -f  *.o *.mod $(executable) core
