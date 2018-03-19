CFLAGS = -Wall
LDLIBS = -lm



.PHONEY : all

all : LSF.png minRosenbrock

LSF.png : LSF.data
	gnuplot LSF.gpi

LSF.data : LSF
	./$< > $@

LSF : LSF.o
	gcc LSF.o -lgsl -lgslcblas -lm -o LSF

LSF.o : LSF.c
	gcc -c LSF.c

minRosenbrock : minRosenbrock.o
	gcc minRosenbrock.o -lgsl -lgslcblas -lm -o minRosenbrock

minRosenbrock.o : minRosenbrock.c
	gcc -c minRosenbrock.c


