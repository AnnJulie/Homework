#ifndef HAVE_QUADRATICSPLINE_H /* This is necessary for multiple includes */
#define HAVE_QUADRATICSPLINE_H

typedef struct {int n; double *x, *y, *b, *c;} qspline;

/* return index of closest point before the value z */
int binary_search(qspline *s, double z);

/* print elements in a vector storing double numbers */
void vector_print(int n, double *v);

/* allocates and builds the quadratic spline */
qspline * qspline_alloc(int n, double *x, double *y, int print);

/* evaluates the prebuilt spline at point z */
double qspline_evaluate(qspline *s, double z);

/* evaluates the derivative of the prebuilt spline at point z */
double qspline_derivative(qspline *s, double z);

/* evaluates the integral of the prebuilt spline from x[0] to z */
double qspline_integral(qspline *s, double z);

/* free memory allocated in qspline_alloc */
void qspline_free(qspline *s);


#endif
