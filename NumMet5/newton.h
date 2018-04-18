#ifndef HAVE_NEWTON_H
#define HAVE_NEWTON_H

#include "QR.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

void newton(void f(gsl_vector* x, gsl_vector* fx), gsl_vector* xstart);
void newton_with_jacobian(void f(gsl_vector* x, gsl_vector* fx, gsl_matrix* J), gsl_vector* xstart);


#endif
