#ifndef HAVE_MINIMIZATION_H
#define HAVE_MINIMIZATION_H

#include "QR.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

/* f: objective function, df: gradient, H: Hessian matrix H*/
void newton(
        double f(gsl_vector* x, gsl_vector* df, gsl_matrix* H),
        gsl_vector* xstart);





#endif
