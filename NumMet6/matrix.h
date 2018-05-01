#ifndef HAVE_MATRIX_H
#define	HAVE_MATRIX_H

#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_linalg.h>

void matrix_generate_random(gsl_matrix* m);
void matrix_print(gsl_matrix* m, char* name);
void matrix_inverse(gsl_matrix* M, gsl_matrix* invM);

#endif
