#ifndef HAVE_INTEGRATOR_H
#define HAVE_INTEGRATOR_H

#include <math.h>
#include <assert.h>
#include <gsl/gsl_vector.h>


double adapt24(double f(double, gsl_vector* calls), double a, double b, double acc, double eps, double f2, double f3, int nrec, gsl_vector* calls);

double adapt(double f(double, gsl_vector* calls), double a, double b, double acc, double eps, gsl_vector* calls);

#endif
