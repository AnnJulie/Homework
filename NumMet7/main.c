#include "ODE.h"
#include "f.h"
#include "matrix.h"
#include <gsl/gsl_matrix.h>

int main()
{
	
	/* setting up parameters for integrating f1 */
	int n=1;		// system is of n equations
	int max = 22;		// max nr of x-values
	double b = 3;	 	// end x-value
	double h = 0.05;	// initial guess for step-size
	double acc = 0.01;	// absolute accuracy goal
	double eps = 0.1;	// relative accuracy goal
	double x0 = 0;		// initial x-value

	gsl_matrix* R = gsl_matrix_calloc(max, n+1); // storing results	


	// evaluating integral
	printf("Evaluating the integral of f(x) = X^2 from 0 to 3:\n");
	double I = ODE_integral(f2, n, x0, b, h, acc, eps, max, R);
	printf("Calculated path\n");
	matrix_print(R, "R");
	printf("I = %g\n", I);

	// free allocated space
	gsl_matrix_free(R);
}

