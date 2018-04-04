CFLAGS += -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS += $$(gsl-config --libs)
CC = gcc

.PHONEY : all clean

all			: main.png
main.png		: main.out.data		;gnuplot main.gpi
main.out.data		: main			;./main < main.in.txt  > $@
main			: main.o linearSpline.o quadraticSpline.o cubicSpline.o
main.o linearSpline.o	: linearSpline.h
main.o quadraticSpline.o: quadraticSpline.h
main.o cubicSpline.o	: cubicSpline.h 


clean:
	-rm -f *.o *.data main
