#ifndef HAVE_ODE_H
#define HAVE_ODE_H

#include <math.h>
#include <gsl/gsl_matrix.h>
#include <stdlib.h>

/* Embedded midpoint-Euler method with error estimate */
void rkstep12(
              void f(int n, double x, double* yx, double* dydx),
              int n,
              double x,
              double* yx,
              double h,
              double* yh,
              double* dy);

/* ODE driver */
int driver(
        void f(int n, double x, double* y, double* dydx),
        int n, double x0, double y0, double b, double h, double acc, double eps, int max,
	gsl_matrix* R);


double ODE_integral(
        void f(int n, double x, double* y, double* dydx),
        int n, double x0, double b, double h, double acc, double eps, int max,
        gsl_matrix* R);

#endif
