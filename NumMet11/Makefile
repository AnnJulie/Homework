CC=gcc-8 
CFLAGS += -Wall -std=gnu1x -O -fopenmp
LDLIBS += -lm -lgomp

.PHONY : out.txt clean

out.txt			: main			;./$< >$@
main			: main.o f.o MC.o
main.o f.o		: f.h
main.o MC.o		: MC.h

clean:
	-rm -f *.o main
