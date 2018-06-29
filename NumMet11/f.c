#include "f.h"
#define PI 3.14159265358979323846


double f1(double *x)
{
	return x[0]*cos(x[0]);
}

double f2(double* x)
{
	return sin(x[0])*sin(x[0]);
}

double f3(double* x)
{
        return cos(x[0])*cos(x[0]);
}
