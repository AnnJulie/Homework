CFLAGS = -Wall -std=gnu99
LDLIBS = -lm

.PHONEY 	: all clean
all		: out.txt		;cat $<		# printing
out.txt 	: main			;./$< > $@	# writing
main 		: main.o komplex.o			# built-in linking
main.o komplex.o: komplex.h				# built-in compiling


clean:
	-rm -f out.txt main main.o komplex.o
