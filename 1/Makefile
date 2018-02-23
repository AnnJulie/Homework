out.txt : hello			# out.txt depends on hello
	./hello > out.txt	# This is how you make out.txt from hello

hello : hello.o			# hello depends on hello.o
	cc hello.o -o hello	# This is how you make hello from hello.o

hello.o : hello.c		# hello.o depends on hello.c
	cc -c hello.c		# This is how you make hello.o from hello.c

.PHONEY:clean
clean:				# "clean" target don't depend on anything
	-rm -f out.txt hello hello.o
