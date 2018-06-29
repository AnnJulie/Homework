#ifndef HAVE_ANN_H
#define HAVE_ANN_H

#include <math.h>
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_multimin.h>

// neural network structure. n is number of hidden neurons.
typedef struct{int n; double (*f)(double); gsl_vector* data;} ann;

// functions
ann* ann_alloc(int n, double (*f)(double));
double ann_feed_forward(ann* network, double x);
void ann_train(ann* network, gsl_vector* xlist, gsl_vector* ylist);
//void ann_train(ann* network, gsl_vector* xlist, gsl_vector* ylist, double delta(ann* network, gsl_vector* xlist, gsl_vector* ylist, const gsl_vector* p, void* params));
void ann_free(ann* network);

double delta(const gsl_vector* p, void* params);

#endif
