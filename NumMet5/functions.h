#ifndef HAVE_FUNCTIONS_H
#define HAVE_FUNCTIONS_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <math.h>

void f1(gsl_vector* x, gsl_vector* fx);
void rosen(gsl_vector* x, gsl_vector* fx);
void himmel(gsl_vector* x, gsl_vector* fx);

void f1_with_J(gsl_vector* x, gsl_vector* fx, gsl_matrix* J);
void rosen_with_J(gsl_vector* x, gsl_vector* fx, gsl_matrix* J);
void himmel_with_J(gsl_vector* x, gsl_vector* fx, gsl_matrix* J);

#endif
