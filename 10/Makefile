CFLAGS = -Wall
LDLIBS = -lm



.PHONEY : all

all : hydrogen.png rosenbrock

hydrogen.png : hydrogen.data
	gnuplot hydrogen.gpi

hydrogen.data : hydrogen
	./$< > $@	

hydrogen : hydrogen.o
	gcc hydrogen.o -lgsl -lgslcblas -lm -o hydrogen

hydrogen.o : hydrogen.c

	gcc -c hydrogen.c

rosenbrock : rosenbrock.o
	gcc rosenbrock.o -lgsl -lgslcblas -lm -o rosenbrock

rosenbrock.o : rosenbrock.c
	gcc -c rosenbrock.c


