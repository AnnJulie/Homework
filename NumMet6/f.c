#include "f.h"


// returning f(x)
double rosenbrock(gsl_vector* var, gsl_vector* df, gsl_matrix* H)
{
        double x = gsl_vector_get(var, 0);
        double y = gsl_vector_get(var, 1);

        // gradient
        gsl_vector_set(df, 0, -2*(1-x) - 400*x*(y-x*x));
        gsl_vector_set(df, 1, 200*(y-x*x));
	// double derivatives
	gsl_matrix_set(H, 0, 0, 2 - 400*y + 1200*x*x);
	gsl_matrix_set(H, 0, 1, -400*x);
	gsl_matrix_set(H, 1, 0, -400*x);
	gsl_matrix_set(H, 1, 1, 200);
	
	return (1-x)*(1-x) + 100*(y-x*x)*(y-x*x);
}


// returning f(x)
double himmelblau(gsl_vector* var, gsl_vector* df, gsl_matrix* H)
{   
	double x = gsl_vector_get(var, 0);
	double y = gsl_vector_get(var, 1);
	// gradient
 	gsl_vector_set(df, 0, 4*x*(x*x+y-11) + 2*(x+y*y-7));
	gsl_vector_set(df, 1, 2*(x*x+y-11) + 4*y*(x+y*y-7));
	// double derivatives
	gsl_matrix_set(H, 0, 0, 12*x*x + 4*y - 42);
        gsl_matrix_set(H, 0, 1, 4*x + 4*y);
        gsl_matrix_set(H, 1, 0, 4*x + 4*y);
        gsl_matrix_set(H, 1, 1, 4*x + 12*y*y - 26);
    
    return (x*x+y-11)*(x*x+y-11) + (x+y*y-7)*(x+y*y-7);
}
