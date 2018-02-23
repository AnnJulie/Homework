#include <stdio.h>
#include "komplex.h"
#define TINY 1e-6


int main()
{
	komplex a = {1,2};
	komplex b = {3,4};

	printf("Testing komplex_add...\n");
	komplex r = komplex_add(a,b);
	komplex R = {4,6};
	komplex_print("a=", a);
	komplex_print("b=", b);
	komplex_print("a+b should = ", R);
	komplex_print("a+b actually = ", r);

	return 0;
}
