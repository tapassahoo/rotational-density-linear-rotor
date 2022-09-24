.SUFFIXES: .o .c .f .cxx

CXX      = g++
CC       = gcc
F77      = gfortran

CFLAGS   = -O3 -DLINUX -DBLAS
CXXFLAGS = -O3 -DLINUX -DBLAS -fopenmp -std=c++11
LDFLAGS  = -lblas -llapack -lm -lgfortran

objects= plgndr.o  
all: run
run: linear_rotor_propagator.o  $(objects) 
	$(CXX) $(CXXFLAGS) -o run linear_rotor_propagator.o $(objects) $(LDFLAGS)
clean:
	rm *.o
.cc.o: 
	$(CXX) -c  $(CXXFLAGS) $<
.c.o:
	$(CC) -c $(LDFLAGS) $<
