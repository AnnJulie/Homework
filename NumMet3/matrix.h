#ifndef HAVE_MATRIX_H
#define	HAVE_MATRIX_H

#include <gsl/gsl_matrix_double.h>
#include <assert.h>

void matrix_generate_random(gsl_matrix* m);
void matrix_print(gsl_matrix* m, char* name);

int matrix_largestElement(gsl_matrix* m, int n);

#endif
