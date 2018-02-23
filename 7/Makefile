CFLAGS = -Wall
LDLIBS = -lm


.PHONEY : all
all : airy.png  linear
	./linear < linear.input.txt > out.io.txt


airy.png : airy.data
	gnuplot airy.gpi

airy.data : airy
	./$< > $@

linear : linear.o
	gcc linear.o -lgsl -lgslcblas -lm -o linear

linear.o : linear.c
	gcc -c linear.c


.PHONEY : clean
clean: 
	rm -f *.o out.io.txt *.data *.png
	
