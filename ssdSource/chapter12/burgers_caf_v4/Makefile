executable=burgers_caf
compile=ftn -h caf
compile=gfortran -fcoarray=single
OPTS=
objects = field.o periodic_2nd_order.o kind_parameters.o 

$(executable): main.o $(objects) Makefile
	$(compile) $(OPTS)  main.o $(objects) -o $(executable) 

main.o: main.F90 $(objects) Makefile
	$(compile) $(OPTS) -c main.F90

periodic_2nd_order.o: periodic_2nd_order.F90 field.o kind_parameters.o Makefile
	$(compile) $(OPTS) -c periodic_2nd_order.F90

field.o: field.F90 kind_parameters.o Makefile
	$(compile) $(OPTS) -c field.F90

kind_parameters.o: kind_parameters.F90 Makefile
	$(compile) $(OPTS) -c kind_parameters.F90

.PHONY : clean
clean :
	-rm -f  *.o *.mod $(executable) core
