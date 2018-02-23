out.txt : math			# out.txt depend on math
	./math > out.txt	# write math to out.txt

math : math.o			# math depend on math.o
	cc math.o -o math	# make math from math.o

math.o : math.c			# math.o depend on math.c
	cc -c math.c		# make math.o from math.c


.PHONEY:clean
clean:
	-rm -f out.txt math math.o
