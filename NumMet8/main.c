#include "integrator.h"
#include "f.h"
#include <math.h>
#include <stdio.h>
#include <gsl/gsl_vector_int.h>

int main()
{	// initializing variables
	gsl_vector* calls = gsl_vector_alloc(1);
	double a = 0;
	double b = 1;
	double acc = 0.001;
	double eps = 0.001;
	double Q;

	// solving intagrals
	printf("\nIntegrating f(x) = sqrt(x) from 0 to 1:\n");
	Q=adapt(f1,a,b,acc,eps, calls);
        printf("Q=%g, calls = %g\n", Q, gsl_vector_get(calls,0));

	printf("\nIntegrating f(x) = 1/sqrt(x) from 0 to 1:\n");
	Q=adapt(f2,a,b,acc,eps, calls);
	printf("Q=%g, calls = %g\n", Q, gsl_vector_get(calls,0));

	printf("\nIntegrating f(x) = ln(x)/sqrt(x) from 0 to 1:\n");
        Q=adapt(f3,a,b,acc,eps, calls);
        printf("Q=%g, calls = %g\n", Q, gsl_vector_get(calls,0));

        printf("\nIntegrating f(x) = 4sqrt(1-(1-x)^2) from 0 to 1:\n");
        Q=adapt(f4,a,b,acc,eps, calls);
        printf("Q=%g, calls = %g\n", Q, gsl_vector_get(calls,0));

	// solving infinite integrals
	double pinf = INFINITY;
	printf("\nIntegrating f(x) = 1/x^2 from 1 to inf:\n");
	Q=adapt(f5,1,pinf,acc,eps, calls);
	printf("Q=%g, calls = %g\n", Q, gsl_vector_get(calls,0));

	return 0;
}

