#include "komplex.h"
#include "stdio.h"

void komplex_print(char* s, komplex z) // prints string s and then komplex z
{
	printf ("%s (%g, %g)\n", s, z.re, z.im);
}


void komplex_set(komplex* z, double x, double y) // z is set to x + iy
{
	(*z).re = x;
	(*z).im = y;
}


komplex komplex_new(double x, double y) // returns x + iy
{
	komplex z = {x,y};
	return z;
}


komplex komplex_add(komplex a, komplex b) // returns a + b
{
	komplex result = {a.re + b.re, a.im + b.im};
	return result;
}


komplex komplex_sub(komplex a, komplex b) // returns a - b
{
	komplex result = {a.re - b.re, a.im - b.im};
	return result;
}


