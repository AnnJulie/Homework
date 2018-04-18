CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

.PHONEY : out.txt clean

out.txt			: main			;./$< > $@
main			: QR.o newton.o functions.o gslFunctions.o matrix.o vector.o
main.o QR.o             : QR.h
main.o newton.o		: newton.h
main.o functions.o      : functions.h
main.o gslFunctions.o  : gslFunctions.h
main.o matrix.o		: matrix.h
main.o vector.o         : vector.h

clean:
	-rm -f *.o main
