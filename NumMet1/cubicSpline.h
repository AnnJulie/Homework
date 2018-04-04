#ifndef HAVE_CUBICSPLINE_H /* This is necessary for multiple includes */
#define HAVE_CUBICSPLINE_H

typedef struct {int n; double *x, *y, *b, *c, *d;} cspline;

/* return index of closest point before the value z */
int binary_search_qubic(cspline *s, double z);

/* allocates and builds the cubic spline */
cspline * cspline_alloc(int n, double *x, double *y);

/* evaluates the prebuilt spline at point z */
double cspline_evaluate(cspline *s, double z);

/* evaluates the derivative of the prebuilt spline at point z */
double cspline_derivative(cspline *s, double z);

/* evaluates the integral of the prebuilt spline from x[0] to z */
double cspline_integral(cspline *s, double z);

/* free memory allocated in qspline_alloc */
void cspline_free(cspline *s);


#endif
