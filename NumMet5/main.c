#include "QR.h"
#include "newton.h"
#include "functions.h"
#include "matrix.h"
#include "vector.h"
#include "gslFunctions.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

int main()
{
	gsl_vector* xstart = gsl_vector_calloc(2);
	
	printf("\n");
	printf("NEWTONS METHOD:\n");

	printf("Finding roots of the linear system.\n");
        printf("Initial guess: x1 = 100, x2 = 1\n");
        gsl_vector_set(xstart, 0, 100);
        gsl_vector_set(xstart, 1, 1);
        newton(f1,xstart);

	printf("Finding roots of the Rosenbrock's valley function.\n");
	printf("Initial guess: x1 = 0.9, x2 = 0.9\n");
	gsl_vector_set(xstart, 0, 0.9);
	gsl_vector_set(xstart, 1, 0.9);	
	newton(rosen,xstart);

	printf("\n");
        printf("Finding roots of the Himmalblau's function.\n");
        printf("Initial guess: x1 = 4, x2 = 3\n");
        gsl_vector_set(xstart, 0, 4);
        gsl_vector_set(xstart, 1, 3); 
        newton(himmel,xstart);


	printf("\n");
        printf("NEWTONS METHOD WITH JACOBIAN:\n");

	printf("Finding roots of the linear system.\n");
        printf("Initial guess: x1 = 100, x2 = 1\n");
        gsl_vector_set(xstart, 0, 100);
        gsl_vector_set(xstart, 1, 1);
        newton_with_jacobian(f1_with_J,xstart);

	printf("Finding roots of the Rosenbrock's valley function.\n");
        printf("Initial guess: x1 = 0.9, x2 = 0.9\n");
        gsl_vector_set(xstart, 0, 0.9);
        gsl_vector_set(xstart, 1, 0.9);
        newton_with_jacobian(rosen_with_J,xstart);

	printf("\n");
        printf("Finding roots of the Himmalblau's function.\n");
        printf("Initial guess: x1 = 4, x2 = 3\n");
        gsl_vector_set(xstart, 0, 4);
        gsl_vector_set(xstart, 1, 3);
        newton_with_jacobian(himmel_with_J,xstart);

	printf("\n");
	printf("GSL procedure for the Rosenbruck function.\n");
        printf("Initial guess: x1 = 0.9, x2 = 0.9\n");
	gsl_multiroot(rosenbrock, rosenbrock_df, rosenbrock_fdf, 0.9, 0.9 );
	
	printf("\n");
	printf("GSL procedure for the Himmelblau function.\n");
        printf("Initial guess: x1 = 4, x2 = 3\n");
        gsl_multiroot(himmelblau, himmelblau_df, himmelblau_fdf, 4, 3 );


	gsl_vector_free(xstart);
}

