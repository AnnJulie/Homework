#ifndef HAVE_GRAMSCHMIDT_H
#define HAVE_GRAMSCHMIDT_H

#include "matrix.h"
#include "vector.h"
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

void qr_gs_inverse(const gsl_matrix* Q, const gsl_matrix* R, gsl_matrix* B);
void qr_gs_decomp(gsl_matrix *A, gsl_matrix* R);
void qr_gs_solve(const gsl_matrix* Q, gsl_matrix* R,const gsl_vector* b, gsl_vector* x);





#endif
