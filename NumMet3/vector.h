#ifndef HAVE_VECTOR_H
#define HAVE_VECTOR_H

#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_matrix_double.h>
#include <assert.h>
#include <math.h>

void vector_get_from_matrix(gsl_vector* v, gsl_matrix* m, size_t col);
void vector_set_to_matrix(gsl_vector* v, gsl_matrix*m, size_t col);

void vector_print(gsl_vector* v, char* name);
double vector_inner_product(gsl_vector* v1, gsl_vector* v2);


#endif
