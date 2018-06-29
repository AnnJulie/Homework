#ifndef HAVE_QSUB_H
#define HAVE_QSUB_H

#include <assert.h>
#include <stdlib.h>
#include <math.h>

typedef struct {int n; double *x, *y, *dy, *c, *d;} qsub_spline;

qsub_spline* qsub_spline_alloc(int n, double *x, double *y, double *dy);
double qsub_spline_eval(qsub_spline* s, double z);
void qsub_spline_free(qsub_spline* s);

#endif
