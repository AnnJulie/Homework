#ifndef HAVE_MC_H
#define HAVE_MC_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void randomx(int dim, double* a, double* b, double* x);

void plainmc(int dim, double* a, double* b, double f(double* x), int N, double* result, double* error);


#endif
