#include "functions.h"


void f1(gsl_vector* x, gsl_vector* fx)
{
	int A = 10000;
	double x0 = gsl_vector_get(x, 0);
	double x1 = gsl_vector_get(x, 1);
	double f1 = A*x0 * x1 - 1;
	double f2 = exp(-1*x0) + exp(-1*x1) - 1 - 1/A;
	gsl_vector_set(fx, 0, f1);
	gsl_vector_set(fx, 1, f2);
}


void rosen(gsl_vector* x, gsl_vector* fx)
{
	double x0 = gsl_vector_get(x, 0);
	double x1 = gsl_vector_get(x, 1);
	double dfdx = -2*(1-x0) - 400*x0*(x1-x0*x0);
	double dfdy = 200*(x1-x0*x0);
	gsl_vector_set(fx, 0, dfdx);
	gsl_vector_set(fx, 1, dfdy);
}


void himmel(gsl_vector* x, gsl_vector* fx)
{
    double x0 = gsl_vector_get(x, 0);
    double x1 = gsl_vector_get(x, 1);
    double dfdx = 4*x0*(x0*x0+x1-11) + 2*(x0+x1*x1-7);
    double dfdy = 2*(x0*x0+x1-11) + 4*x1*(x0+x1*x1-7);
    gsl_vector_set(fx, 0, dfdx);
    gsl_vector_set(fx, 1, dfdy);
}



void f1_with_J(gsl_vector* x, gsl_vector* fx, gsl_matrix* J)
{
    int A = 10000;
    double x0 = gsl_vector_get(x, 0);
    double x1 = gsl_vector_get(x, 1);
    double f1 = A*x0 * x1 - 1;
    double f2 = exp(-1*x0) + exp(-1*x1) - 1 - 1/A;
    gsl_vector_set(fx, 0, f1);
    gsl_vector_set(fx, 1, f2);
    gsl_matrix_set(J, 0, 0, A*x1);
    gsl_matrix_set(J, 0, 1, A*x0);
    gsl_matrix_set(J, 1, 0, -x0);
    gsl_matrix_set(J, 1, 1, -x1);
}


void rosen_with_J(gsl_vector* x, gsl_vector* fx, gsl_matrix* J)
{
        double x0 = gsl_vector_get(x, 0);
        double x1 = gsl_vector_get(x, 1);
        double dfdx = -2*(1-x0) - 400*x0*(x1-x0*x0);
        double dfdy = 200*(x1-x0*x0);
        gsl_vector_set(fx, 0, dfdx);
        gsl_vector_set(fx, 1, dfdy);
	gsl_matrix_set(J, 0, 0, 2 - 400*x1 + 1200*x0*x0);
	gsl_matrix_set(J, 0, 1, -400*x0);
	gsl_matrix_set(J, 1, 0, -400*x0);
	gsl_matrix_set(J, 1, 1, 200);
}

void himmel_with_J(gsl_vector* x, gsl_vector* fx, gsl_matrix* J)
{   
	double x0 = gsl_vector_get(x, 0);
	double x1 = gsl_vector_get(x, 1);
	double dfdx = 4*x0*(x0*x0+x1-11) + 2*(x0+x1*x1-7);
	double dfdy = 2*(x0*x0+x1-11) + 4*x1*(x0+x1*x1-7);
 	gsl_vector_set(fx, 0, dfdx);
	gsl_vector_set(fx, 1, dfdy);
	gsl_matrix_set(J, 0, 0, 12*x0*x0 + 4*x1 - 42);
        gsl_matrix_set(J, 0, 1, 4*x0 + 4*x1);
        gsl_matrix_set(J, 1, 0, 4*x0 + 4*x1);
        gsl_matrix_set(J, 1, 1, 4*x0 + 12*x1*x1 - 26);
}
