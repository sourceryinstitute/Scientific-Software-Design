#The following compilers have bugs that prevent compiling this code correctly
#(Bug reports have been submitted to the compiler team):
compile=gfortran -c
#The following compilers compile this code correctly:
compile=xlf2003 -c
compile=nagfor -f2003

universal.o: universal.F90 hermetic.o ref_counter.o Makefile
	$(compile) -c universal.F90

ref_counter.o: ref_counter.F90 hermetic.o Makefile
	$(compile) -c ref_counter.F90

hermetic.o: hermetic.F90 Makefile
	$(compile) -c hermetic.F90

.PHONY : clean
clean :
	-rm -f  *.o *.mod $(executable)
