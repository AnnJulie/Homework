#ifndef	HAVE_JACOBIDIAGONALIZATION_H
#define HAVE_JACOBIDIAGONALIZATION_H

#include <math.h>
#include "matrix.h"
#include "vector.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <time.h>

int jacobiCyclic(gsl_matrix* A, gsl_vector* e, gsl_matrix* V);
int lowestEigenvaluesCyclic(gsl_matrix* A, gsl_vector* e, gsl_matrix* V);
int jacobiValueByValue(gsl_matrix* A, gsl_vector* e, gsl_matrix* V);

double cyclicDiagTime(size_t n);
double valueDiagTime(size_t n);

void testCyclicDiagonalization();


#endif
