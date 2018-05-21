#include "f.h"
#include <stdio.h>

double f1(double x, gsl_vector* calls)
{
	gsl_vector_set(calls, 0, gsl_vector_get(calls, 0)+1);
	return sqrt(x);
}

double f2(double x, gsl_vector* calls)
{
        gsl_vector_set(calls, 0, gsl_vector_get(calls, 0)+1);
        return 1/sqrt(x);
}

double f3(double x, gsl_vector* calls)
{
        gsl_vector_set(calls, 0, gsl_vector_get(calls, 0)+1);
        return log(x)/sqrt(x);
}

double f4(double x, gsl_vector* calls)
{
        gsl_vector_set(calls, 0, gsl_vector_get(calls, 0)+1);
        return 4*sqrt(1-(1-x)*(1-x));
}

double f5(double x, gsl_vector* calls)
{
        gsl_vector_set(calls, 0, gsl_vector_get(calls, 0)+1);
        return 1/(x*x);
}
