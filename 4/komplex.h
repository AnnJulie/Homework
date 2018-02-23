#ifndef HAVE_KOMPLEX_H /* This is necessary for multiple includes */
#define HAVE_KOMPLEX_H

typedef struct {double re, im;} komplex; // define a structure of complex nrs

void komplex_print	(char* s, komplex z); // prints string s and then komplex z
void komplex_set	(komplex* z, double x, double y); // z is set to x + iy
komplex komplex_new	(double x, double y); // returns x + iy
komplex komplex_add	(komplex a, komplex b); // returns a + b
komplex komplex_sub	(komplex a, komplex b); // returns a - b

#endif
