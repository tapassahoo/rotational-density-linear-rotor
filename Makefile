.SUFFIXES: .o .c .f .cxx

options= -O3 -DLINUX -DBLAS
CC=g++ #usr/bin/g++ #/opt/openmpi/bin/mpic++
cc=gcc #/usr/bin/gcc #/opt/openmpi/bin/mpicc
f77=gfortran #/opt/openmpi/bin/mpif90

objects= plgndr.o  
all: run
run: linear_rotor_propagator.o  $(objects) 
	$(CC) $(options) -o run linear_rotor_propagator.o $(objects)  -L/usr/local/lib -lgfortran
clean:
	rm *.o
.cxx.o: 
	$(CC) -c  $(options) $<
.c.o:
	$(cc) -c $(options) $<
