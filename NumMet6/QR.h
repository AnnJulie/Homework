#ifndef HAVE_QR_H
#define HAVE_QR_H

#include "matrix.h"
#include "vector.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>


void backSubstitute(gsl_matrix* A, gsl_vector* x, gsl_vector* b);
void gramschmidt(gsl_matrix* A, gsl_matrix* R);
void QR_solve(gsl_matrix* A, gsl_vector* x, gsl_vector* b);
void QR_generate_example(gsl_matrix* A, gsl_vector* b);


#endif
