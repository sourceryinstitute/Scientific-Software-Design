include_directories(../../chapter06/common_cpp/)
add_executable(timed_lorenz_cpp explicit_euler.cpp integrand.cpp lorenz.cpp main.cpp runge_kutta_2nd.cpp timed_lorenz.cpp)
