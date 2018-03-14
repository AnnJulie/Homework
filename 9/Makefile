CFLAGS = -Wall
LDLIBS = -lm



.PHONEY : all

all : integral varMethod.png

varMethod.png : varMethod.data
	gnuplot varMethod.gpi

varMethod.data : varMethod
	./$< > $@	

varMethod : varMethod.o
	gcc varMethod.o -lgsl -lgslcblas -lm -o varMethod

varMethod.o : varMethod.c
	gcc -c varMethod.c

integral : integral.o
	gcc integral.o -lgsl -lgslcblas -lm -o integral

integral.o : integral.c
	gcc -c integral.c

