CFLAGS += -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS += $$(gsl-config --libs)

.PHONY : all clean


all			: main.png
main.png		: main.out.txt		;gnuplot main.gpi
main.out.txt		: main			;./main< main.in.txt >$@
main			: main.o qsub.o
main.o qsub.o           : qsub.h


clean:
	-rm -f *.o main
