#include <stdio.h>
#include <math.h>
#include <gsl/gsl_multimin.h>

/* return the function evaluated in (x,y) */
double rosenbrock(const gsl_vector *v, void *params)
{
	double x,y;
	x = gsl_vector_get(v, 0);
	y = gsl_vector_get(v, 1);

	return (1-x)*(1-x) + 100*(y-x*x)*(y-x*x);
}

/* stores the partial derivatives in df */
void rosenbrock_df(const gsl_vector *v, void *params, gsl_vector *df)
{
	double x,y;
	x = gsl_vector_get (v, 0);
	y = gsl_vector_get (v, 1);

	gsl_vector_set (df, 0, -2*(1-x)-400*x*(y-x*x));
	gsl_vector_set (df, 1, 200*(y-x*x));
}

/* optimisation of the two above procedures by computing them at the same time */
void rosenbrock_fdf (const gsl_vector *v, void *params, double *f, gsl_vector *df)
{
	*f = rosenbrock(v, params);
	rosenbrock_df(v, params, df);
}


int main(void)
{
	size_t iter = 0;
	int status;

	const gsl_multimin_fdfminimizer_type *T;
	gsl_multimin_fdfminimizer *s;

	gsl_vector *x;
	gsl_multimin_function_fdf my_func;

	my_func.n = 2;
	my_func.f = rosenbrock;
	my_func.df = rosenbrock_df;
	my_func.fdf = rosenbrock_fdf;
	my_func.params = NULL;	

	/* starting point x,y = (0,0) */
	x = gsl_vector_alloc (2);
	gsl_vector_set (x, 0, 0);
	gsl_vector_set (x, 1, 0);

	T = gsl_multimin_fdfminimizer_conjugate_fr;
	s = gsl_multimin_fdfminimizer_alloc (T, 2);

	gsl_multimin_fdfminimizer_set (s, &my_func, x, 0.01, 1e-4);

	do{
		iter++;
		status = gsl_multimin_fdfminimizer_iterate (s);
		if(status){break;}

		status = gsl_multimin_test_gradient (s->gradient, 1e-3);
		if(status == GSL_SUCCESS){
			printf("\nminimum found! \n");
			printf("iterations: %5zu\n", iter);
			printf("(x,y) = (%.5f %.5f)\n\n",gsl_vector_get(s->x, 0), gsl_vector_get(s->x, 1));
		}
	}
	while (status == GSL_CONTINUE && iter < 100);

	gsl_multimin_fdfminimizer_free (s);
	gsl_vector_free (x);

	return 0;
}

