CFLAGS = -Wall -std=gnu99
LDLIBS = -lm
CC = gcc

.PHONEY: all
all: out.txt
	cat out.txt

out.txt : main
	./main > out.txt

main : main.o equal.o
	gcc main.o equal.o -o main

main.o : main.c
	gcc -c main.c -o main.o

equal.o : equal.c
	gcc -c equal.c -o equal.o



.PHONEY:clean
clean:
	-rm -f out.txt main main.o equal.o
