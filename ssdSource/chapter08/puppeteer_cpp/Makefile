executable=chaos_cpp
compile=/gsa/tlbgsa/projects/x/xlcmpbld/run/vacpp/dev/aix/daily/latest/bin/xlC -qrtti=all
compile=g++
INCLUDE=-I../../chapter06/common_cpp/

$(executable): main.o air.o cloud.o ground.o atmosphere.o integrand.o mat.o
	$(compile) main.o air.o cloud.o ground.o atmosphere.o integrand.o mat.o -o $(executable) 

main.o: main.cpp 
	$(compile) $(INCLUDE) -c main.cpp 

integrand.o: integrand.cpp integrand.h
	$(compile) $(INCLUDE) -c integrand.cpp 

air.o: air.cpp air.h
	$(compile) $(INCLUDE) -c air.cpp 

cloud.o: cloud.cpp cloud.h
	$(compile) $(INCLUDE) -c cloud.cpp 

ground.o: ground.cpp ground.h
	$(compile) $(INCLUDE) -c ground.cpp 

atmosphere.o: atmosphere.cpp atmosphere.h
	$(compile) $(INCLUDE) -c atmosphere.cpp 

mat.o: mat.cpp mat.h
	$(compile) $(INCLUDE) -c mat.cpp 

.PHONY : clean
clean :
	-rm -f  *.o *.mod $(executable)
