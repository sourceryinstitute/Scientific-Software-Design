executable=burgers_cpp
INCLUDE=-I../../chapter06/common_cpp/ -I../../chapter08/puppeteer_cpp/ -I../../appendixA/
compile=/gsa/tlbgsa/projects/x/xlcmpbld/run/vacpp/dev/aix/daily/latest/bin/xlC $(INCLUDE) -qrtti=all
compile=g++
objects = periodic_6th_factory.o periodic_6th_order.o initializer.o mat.o gaussian_elimination.o

$(executable): main.o $(objects)
	$(compile)  main.o $(objects) -o $(executable)

main.o: main.cpp $(objects)
	$(compile) $(INCLUDE) -c main.cpp

periodic_6th_factory.o: periodic_6th_factory.h periodic_6th_factory.cpp periodic_6th_order.o field_factory.h
	$(compile) $(INCLUDE) -c periodic_6th_factory.cpp

periodic_6th_order.o: periodic_6th_order.h periodic_6th_order.cpp 
	$(compile) $(INCLUDE) -c periodic_6th_order.cpp

initializer.o: initializer.h initializer.cpp
	$(compile) $(INCLUDE) -c initializer.cpp

mat.o: ../../chapter08/puppeteer_cpp/mat.h ../../chapter08/puppeteer_cpp/mat.cpp
	$(compile) $(INCLUDE) -c ../../chapter08/puppeteer_cpp/mat.cpp

gaussian_elimination.o: ../../appendixA/gaussian_elimination.h ../../appendixA/gaussian_elimination.cpp
	$(compile) $(INCLUDE) -c ../../appendixA/gaussian_elimination.cpp

.PHONY : clean
clean :
	-rm -f  *.o $(executable) core
