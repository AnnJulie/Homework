#include "QR.h"
#include "matrix.h"
#include "vector.h"
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_blas.h>

int main()
{
	int fns = 3;
	int n = 9;
	double x[] = {-2.1, -1.5, -1.03, -0.42, 0.02, 0.58, 1.06, 1.45, 2.03};
	double y[] = {2.1, 0.73, -0.25, -0.65, 0.04, 0.85, 2.06, 3.65, 6.15};
	double dy[]= {1.24, 0.894, 1.31, 0.898, 0.71, 0.598, 0.535, 0.768, 0.978};
	fit_data* data = fit_data_alloc(n, x, y, dy);

	gsl_matrix* S = gsl_matrix_calloc(fns,fns);
	gsl_vector* c = gsl_vector_calloc(fns);
	lsfit(fns, data, c, S);
	
	// write fit-function to file
	double c0 = gsl_vector_get(c,0);
	double c1 = gsl_vector_get(c,1);
	double c2 = gsl_vector_get(c,2);
	double dc0 = gsl_matrix_get(S,0,0);
	double dc1 = gsl_matrix_get(S,1,1);
	double dc2 = gsl_matrix_get(S,2,2);

	for (int i=0; i<n; i++){printf("%g %g %g\n", x[i], y[i], dy[i]);}
	printf("\n\n");
	for (double i=x[0]; i<x[n-1]; i+=0.1){
		printf("%g %g %g %g\n", 
			i, 
			c0*funs(0,i) + c1*funs(1,i) + c2*funs(2,i), 
			(c0+dc0)*funs(0,i) + (c1+dc1)*funs(1,i) + (c2+dc2)*funs(2,i),
			(c0-dc0)*funs(0,i) + (c1-dc1)*funs(1,i) + (c2-dc2)*funs(2,i) );
	}	

	printf("\n\n");
	vector_print(c, "c");
	matrix_print(S, "S");

	fit_data_free(data);
	gsl_matrix_free(S);
	gsl_vector_free(c);
}

