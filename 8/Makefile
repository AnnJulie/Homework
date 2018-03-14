CFLAGS = -Wall
LDLIBS = -lm



.PHONEY : orbit.png

orbit.png : orbit.data
	gnuplot orbit.gpi

orbit.data : orbit Makefile
	./orbit -p -0.7 -e 0.02 > $@

orbit : orbit.o
	gcc orbit.o -lgsl -lgslcblas -lm -o orbit

orbit.o : orbit.c
	gcc -c orbit.c

diff.png : diff.data analytical.data
	gnuplot diff.gpi

diff.data : diff
	./$< > $@

diff : diff.o
	gcc diff.o -lgsl -lgslcblas -lm -o diff

diff.o : diff.c
	gcc -c diff.c

analytical.data : analytical
	./$< > $@
