CFLAGS = -Wall
LDLIBS = -lm

.PHONEY : error.pdf

error.pdf : error.data
	gnuplot error.gpi

error.data : error
	./error -3 3 0.2 > error.data

error : error.o
	gcc error.o -lgsl -lgslcblas -lm -o error

error.o : error.c
	gcc -c error.c

.PHONEY : clean
clean: 
	-rm -f error.o error.data error
