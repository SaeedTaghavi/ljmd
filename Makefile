################
# Compiler settings and flags
CC=g++
#OPENMP=-fopenmp
OPENMP=
CFLAGS= -O0 -Wall -W -g $(OPENMP)
#CFLAGS=-O3 -march=native -ffast-math -fexpensive-optimizations -ftree-vectorize $(OPENMP)
################

# $@ is the name of the file being generated (ljmd)
# $^ are the source files (mdsys.cpp, ljmd.cpp)
ljmd: mdsys.cpp ljmd.cpp
	$(CC) -o $@ $(CFLAGS) $^ -lm

clean:
	rm -f ljmd *.dat *.xyz *.o *.mod perf.data

bench: bench1 bench2 bench3

bench1: ljmd
	time ./ljmd argon_108.inp

bench2:
	time ./ljmd argon_2916.inp

bench3:
	time ./ljmd argon_78732.inp

.PHONY: clean bench bench1 bench2 bench3
