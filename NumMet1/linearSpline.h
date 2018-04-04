#ifndef HAVE_LINEARSPLINE_H /* This is necessary for multiple includes */
#define HAVE_LINEARSPLINE_H


double linterp(int n, double *x, double *y, double z);

double linterp_integ(int n, double *x, double *y, double z);

#endif
