CC = gcc-8
CFLAGS += -Wall -std=gnu1x $$(gsl-config --cflags) -Ofast
LDLIBS += $$(gsl-config --libs)

.PHONY : all clean

all			: annFit.png
annFit.png		: annFit.data		; gnuplot annFit.gpi
annFit.data		: main			;./$< >$@
main			: main.o ann.o
main.o ann.o		: ann.h

clean:
	-rm -f *.o main
