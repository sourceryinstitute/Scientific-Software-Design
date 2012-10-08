include_directories(../../chapter06/common_cpp/)
include_directories(../../appendixA/)
include_directories(../../chapter08/puppeteer_cpp/)

add_executable(burgers_cpp main.cpp periodic_6th_factory.cpp periodic_6th_order.cpp initializer.cpp ../../appendixA/gaussian_elimination.cpp ../../chapter08/puppeteer_cpp/mat.cpp)
