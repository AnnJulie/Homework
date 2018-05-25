#include "f.h"
#define PI 3.14159265358979323846

double f1(double* x)
{
	return x[0]*x[1];
}

double f2(double* x)
{
        return 1/(PI*PI*PI*(1-cos(x[0])*cos(x[1])*cos(x[2])));
}
