executable=chaos_cpp
compile=/gsa/tlbgsa/projects/x/xlcmpbld/run/vacpp/dev/aix/daily/latest/bin/xlC -qrtti=all
compile=g++
INCLUDE=-I../../chapter06/common_cpp/

$(executable): main.o integrand.o lorenz.o timed_lorenz.o runge_kutta_2nd.o explicit_euler.o
	$(compile) main.o integrand.o lorenz.o timed_lorenz.o runge_kutta_2nd.o explicit_euler.o -o $(executable) 

main.o: main.cpp 
	$(compile) $(INCLUDE) -c main.cpp 

runge_kutta_2nd.o: runge_kutta_2nd.cpp runge_kutta_2nd.h
	$(compile) $(INCLUDE) -c runge_kutta_2nd.cpp 

explicit_euler.o: explicit_euler.cpp explicit_euler.h
	$(compile) $(INCLUDE) -c explicit_euler.cpp 

timed_lorenz.o: timed_lorenz.cpp timed_lorenz.h 
	$(compile) $(INCLUDE) -c timed_lorenz.cpp

lorenz.o: lorenz.cpp lorenz.h 
	$(compile) $(INCLUDE) -c lorenz.cpp

integrand.o: integrand.cpp integrand.h
	$(compile) $(INCLUDE) -c integrand.cpp 

.PHONY : clean
clean :
	-rm -f  *.o *.mod $(executable)
