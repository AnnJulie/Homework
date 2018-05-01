#ifndef HAVE_F_H
#define HAVE_F_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <math.h>

double rosenbrock(gsl_vector* x, gsl_vector* df, gsl_matrix* H);
double himmelblau(gsl_vector* x, gsl_vector* df, gsl_matrix* H);

#endif
