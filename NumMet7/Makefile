CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

.PHONEY : main.out clean

main.out		: main			;./$< >$@
main			: ODE.o f.o matrix.o
main.o ODE.o            : ODE.h
main.o f.o              : f.h
main.o matrix.o         : matrix.h

clean:
	-rm -f *.o main
