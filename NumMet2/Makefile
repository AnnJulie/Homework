CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

.PHONEY : all clean

all			: main			;./main
main			: gramSchmidt.o matrix.o vector.o
main.o gramSchmidt.o	: gramSchmidt.h
main.o matrix.o		: matrix.h
main.o vector.o         : vector.h

clean:
	-rm -f *.o main
