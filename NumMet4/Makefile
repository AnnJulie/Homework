CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

.PHONEY : all clean

all			: lsfit.png
lsfit.png		: lsfit.data		;gnuplot lsfit.gpi
lsfit.data		: main			;./$< > $@
main			: QR.o matrix.o vector.o
main.o QR.o		: QR.h
main.o matrix.o		: matrix.h
main.o vector.o         : vector.h

clean:
	-rm -f *.o *.data main
