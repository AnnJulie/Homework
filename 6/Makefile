CFLAGS = -Wall -std=gnu99
LDLIBS = -lm
FORMAT = pdf

.PHONEY : all
all : main-stdio main-cmdline
	gcc main-stdio.o -o main-stdio
	gcc main-cmdline.o -o main.cmdline


main-stdio.o : main-stdio.c
	gcc -c main-stdio.c

main-cmdline.o : main-cmdline.c
	gcc -c main-cmdline.c


.PHONEY : clean
clean:
	-rm -f main-stdio main-stdio.o main-cmdline main-cmdline.o
