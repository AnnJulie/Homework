CFLAGS = -Wall -std=c99
LDLIBS = -lm
CC = clang
CC = gcc


.PHONEY		: all clean

all : out.txt
	cat out.txt		# printing

out.txt	: main
	./main > out.txt	# writing

main: main.o nvector.o		# built-in linking

main.o nvector.o: nvector.h	# build-in compiling

clean:
	-rm -f out.txt main main.o nvector.o	
