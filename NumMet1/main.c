#include "linearSpline.h"
#include "quadraticSpline.h"
#include "cubicSpline.h"
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_spline.h>

int main(){
        int n = 24;
        double z = 2.5;
        double x[n], y[n], e1, e2;
        int i = 0;
        /* reading input table, assign to vectors, write table to out-file */
        while (scanf("%lg", &e1) != EOF && scanf("%lg", &e2) != EOF && i<n){
                x[i] = e1;
                y[i] = e2;
                printf("%g %g\n", e1, e2);
                i++;
        }

        /* linear spline */
        printf("\n\n"); /* Gnuplot index 1 */
        printf("%g %g\n", z, linterp(n, x, y, z));
	fprintf(stderr, "\nLinear spline: f(z) = %g\n", linterp(n, x, y, z));
        fprintf(stderr, "Integral from x[0] to z: I = %g\n", linterp_integ(n, x,y,z));

	/* quadratic spline */
	qspline* s = qspline_alloc(n, x, y, 0);
	fprintf(stderr, "\nQuadratic spline: f(z) = %g\n", qspline_evaluate(s, z));
	fprintf(stderr, "Derivative: f'(z) = %g\n", qspline_derivative(s, z));
	fprintf(stderr, "Integral from x[0] to z: I = %g\n", qspline_integral(s, z));
	qspline_free(s);

	/* Test ing the code with some easy tables */
	int test_n = 5;	
	double test_X[] = {1,2,3,4,5};
	double test_Y1[] = {1,1,1,1,1};
	double test_Y2[] = {1,2,3,4,5};
	double test_Y3[] = {1,4,9,16,25};

	printf("\n\nTest 1:\n");
	qspline* s_test1 = qspline_alloc(test_n, test_X, test_Y1, 1);

	printf("\n\nTest 2:\n");
        qspline* s_test2 = qspline_alloc(test_n, test_X, test_Y2, 1);
	
	printf("\n\nTest 3:\n");
        qspline* s_test3 = qspline_alloc(test_n, test_X, test_Y3, 1);

	qspline_free(s_test1);
	qspline_free(s_test2);
	qspline_free(s_test3);


	/* Cubic Spline */
	cspline* cs = cspline_alloc(n, x, y);
	fprintf(stderr, "\nCubic spline: f(z) = %g\n", cspline_evaluate(cs, z));
	fprintf(stderr, "Derivative: f'(z) = %g\n", cspline_derivative(cs, z));
	fprintf(stderr, "Integral from x[0] to z: I = %g\n", cspline_integral(cs, z));        
	cspline_free(cs);

	/* comparing with GSL cubic spline functions */
	const gsl_spline* gs = gsl_spline_alloc(gsl_interp_cspline, 24);
	int init = gsl_spline_init(gs, &x[0], &y[0], 24);

	gsl_interp_accel* a = gsl_interp_accel_alloc();
	fprintf(stderr, "\nGSL spline: f(z) = %g\n", gsl_spline_eval(gs, z, a));
	fprintf(stderr, "Derivative: f'(z) = %g\n", gsl_spline_eval_deriv(gs, z, a));
	fprintf(stderr, "Integral from x[0] to z: I = %g\n", gsl_spline_eval_integ(gs, x[0], z, a));
}


