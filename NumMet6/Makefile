CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

.PHONEY : main.out clean

main.out		: main			;./$< > $@
main			: f.o QR.o matrix.o vector.o minimization.o
main.o f.o              : f.h
main.o QR.o             : QR.h
main.o matrix.o		: matrix.h
main.o vector.o         : vector.h
main.o minimization.o   : minimization.h
clean:
	-rm -f *.o main
