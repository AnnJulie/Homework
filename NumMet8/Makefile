CFLAGS = -Wall -std=c99 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

.PHONEY : main.out clean

main.out		: main			;./$< >$@
main			: integrator.o f.o
main.o integrator.o     : integrator.h
main.o f.o              : f.h


clean:
	-rm -f *.o main
