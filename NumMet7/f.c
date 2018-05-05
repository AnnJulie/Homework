#include "f.h"


void f1(int n, double x, double* yx, double* dydx)
{
	dydx[0] = 2-exp(-4*x) -2*yx[0];
}


void f2(int n, double x, double* yx, double* dydx)
{
	dydx[0] = x*x;
}
