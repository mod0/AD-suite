CXX=g++
CXXFLAGS=-O3 -std=c++11

all : stream_vel check_gradient check_hessian

stream_vel : stream_vel.cpp stream_vel.hpp
	$(CXX) $(CXXFLAGS) -o $@ stream_vel.cpp

check_gradient : stream_vel.cpp stream_vel.hpp
	$(CXX) $(CXXFLAGS) -D CHECK_GRADIENT -o $@ stream_vel.cpp

check_hessian : stream_vel.cpp stream_vel.hpp
	$(CXX) $(CXXFLAGS) -D CHECK_HESSIAN -o $@ stream_vel.cpp

clean:
	rm -f stream_vel check_gradient check_hessian

.PHONY: clean
