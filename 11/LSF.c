#include <stdio.h>
#include <math.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>

struct exp_dat {int n; double *t, *y, *e;};

double master (const gsl_vector *v, void *params)
{
	double A = gsl_vector_get(v, 0);
	double T = gsl_vector_get(v, 1);
	double B = gsl_vector_get(v, 2);
	struct exp_dat *p = (struct exp_dat*) params;
	
	int n = p->n;
	double *t = p->t;
	double *y = p->y;
	double *e = p->e;

	double sum = 0;
	#define f(t) A*exp(-(t)/T) + B
	for (int i=0; i<n; i++){
		sum += pow(  (f(t[i]) - y[i]) / e[i], 2);
	}
	return sum;
}


int main(void)
{
	/* Initialize data */
	int n = 10;
        double t[] = {0.47,1.41,2.36,3.30,4.24,5.18,6.13,7.07,8.01,8.95};
        double y[] = {5.49,4.08,3.54,2.61,2.09,1.91,1.55,1.47,1.45,1.25};
        double e[] = {0.26,0.12,0.27,0.10,0.15,0.11,0.13,0.07,0.15,0.09};

	struct exp_dat p = {.n=n, .t=t, .y=y, .e=e};

	/* Setting up the environment */
	const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
	gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc(T,3);
	gsl_vector *ss, *x;
	gsl_multimin_function min_func;

        min_func.n = 3;
        min_func.f = master;
        min_func.params = (struct exp_dat*) &p;

	size_t iter = 0;
	int status;
	double size;

	/* Starting point */
	x = gsl_vector_alloc (3);
	gsl_vector_set (x, 0, 6.0);
	gsl_vector_set (x, 1, 1.0);
	gsl_vector_set (x, 2, 0.1);	

	/* Set initial step sizes */
	ss = gsl_vector_alloc (3);
	gsl_vector_set (ss, 0, 0.1);
	gsl_vector_set (ss, 1, 0.1);
	gsl_vector_set (ss, 2, 0.1);


	/* iterate */
	gsl_multimin_fminimizer_set (s, &min_func, x, ss);
	do{
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);
		if(status) {break;}

		size = gsl_multimin_fminimizer_size (s);
		status = gsl_multimin_test_size (size, 1e-2);
	
		if (status == GSL_SUCCESS)
		{
			fprintf(stderr,"Fitting Converged!\n\n");
		}		
	}while (status == GSL_CONTINUE && iter < 100);

	double Af = gsl_vector_get (s->x, 0);
	double Tf = gsl_vector_get (s->x, 1);
	double Bf = gsl_vector_get (s->x, 2);
	
	/* writing parameters for plotting */
	for (int i=0; i<n; i++){
		printf("%g %g\n", t[i], y[i]);
	}
	printf("\n\n");
	
	double tt;
	#define ff(t) Af*exp(-(t)/Tf) + Bf
	for (tt=0; tt<10; tt+=0.1){
		printf("%g %g\n", tt, ff(tt));
	}

	/* free allocated space */
	gsl_vector_free(x);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free (s);

	return status;
}	



