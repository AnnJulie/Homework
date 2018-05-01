#include "minimization.h"
#include "f.h"
#include "QR.h"
#include "matrix.h"
#include "vector.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>


int main()
{
	gsl_vector* xstart = gsl_vector_calloc(2);	

// -----ROSENBROCK---------------//
	printf("\n");
        printf("Minimize the Rosenbrock's valley function.\n");
        printf("Initial guess: x = 0.5, y = 1.5\n");
        gsl_vector_set(xstart, 0, 0.5);
        gsl_vector_set(xstart, 1, 1.5);
        newton(rosenbrock, xstart);

// -----Himmelblau---------------//
	printf("\n");
	printf("Minimize the Himmelblau's function.\n");
	printf("Initial guess: x = 3.5, y = 1.7\n");
	gsl_vector_set(xstart, 0, 3.5);
	gsl_vector_set(xstart, 1, 1.7);	
	newton(himmelblau, xstart);

//------Fitting-problem----------//
double t[] = {0.23,1.29,2.35,3.41,4.47,5.53,6.59,7.65,8.71,9.77};
double y[] = {4.64,3.38,3.01,2.55,2.29,1.67,1.59,1.69,1.38,1.46};
double e[] = {0.42,0.37,0.34,0.31,0.29,0.27,0.26,0.25,0.24,0.24};
int N = sizeof(t)/sizeof(t[0]);



// free allocated space
gsl_vector_free(xstart);
}

