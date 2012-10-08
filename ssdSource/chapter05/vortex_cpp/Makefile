executable=quantumvortex
objects=vortex.o
compile=/gsa/tlbgsa/projects/x/xlcmpbld/run/vacpp/111/aix/daily/latest/bin/xlC -qrtti=all
compile=/opt/local/bin/g++-mp-4.6
compile=g++
INCLUDE=-I../../chapter06/common_cpp/

$(executable): main.o $(objects)
	$(compile) main.o $(objects) -o $(executable) 

main.o: main.cpp vortex.o
	$(compile) $(INCLUDE) -c main.cpp

vortex.o: vortex.cpp
	$(compile) $(INCLUDE) -c vortex.cpp 

.PHONY : clean
clean :
	-rm -f  *.o $(executable) a.out core
