executable=chaos
compile=/gsa/tlbgsa/projects/x/xlcmpbld/run/vacpp/dev/aix/daily/latest/bin/xlC -qrtti=all
compile=g++
INCLUDE=-I../common_cpp/

$(executable): main.o integrand.o lorenz.o 
	$(compile) main.o integrand.o lorenz.o -o $(executable) 

main.o: main.cpp lorenz.o 
	$(compile) $(INCLUDE) -c main.cpp 

lorenz.o: lorenz.cpp lorenz.h integrand.o
	$(compile) $(INCLUDE) -c lorenz.cpp

integrand.o: integrand.cpp integrand.h
	$(compile) $(INCLUDE) -c integrand.cpp 

.PHONY : clean
clean :
	-rm -f  *.o *.mod $(executable)
