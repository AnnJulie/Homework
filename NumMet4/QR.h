#ifndef HAVE_QR_H
#define HAVE_QR_H

#include "matrix.h"
#include "vector.h"
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

typedef struct {int n; double *x, *y, *dy;} fit_data;

fit_data* fit_data_alloc(int n, double* x, double* y, double* dy);
void fit_data_free(fit_data* data);
double funs(int i, double x);
void lsfit(int fs, fit_data* data, gsl_vector* c, gsl_matrix* S);
void backSubstitute(gsl_matrix* A, gsl_vector* x, gsl_vector* b);
void gramschmidt(gsl_matrix* Q, gsl_matrix* R);
void QR_generate_example(gsl_matrix* A, gsl_vector* b);

#endif
