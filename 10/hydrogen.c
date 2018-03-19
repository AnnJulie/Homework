#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multiroots.h>


int radial_f(double r, const double f[], double dfdr[], void *params)
{
	(void)(r); /* avoid unused parameters */

	double e = *(double *) params;

	dfdr[0] = f[1];
	dfdr[1] = -2*(1/r+e)*f[0];
	return GSL_SUCCESS;
}

double Fe(double e, double r)
{
	assert(r>=0.0);
	const double rmin = 1e-3;
	if(r<rmin){return r-r*r;}

	gsl_odeiv2_system sys ={radial_f, NULL, 2, &e};

	gsl_odeiv2_driver *driver =
		gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, 1e-3, 1e-6, 1e-6);

	double t = rmin;
	double f[2] = {t-t*t, 1-2*t};

	int status = gsl_odeiv2_driver_apply(driver, &t, r, f);
	if(status != GSL_SUCCESS) {fprintf(stderr, "shit! error in fe: %d\n", status);}

	gsl_odeiv2_driver_free(driver);

	return f[0];
}


int master(const gsl_vector *ev, void *params, gsl_vector *f){
	double e = gsl_vector_get(ev, 0);
	assert(e<0);
	double rmax = *(double*)params;
	double fval = Fe(e, rmax);
	gsl_vector_set(f, 0, fval);
	return GSL_SUCCESS;
}



int main()
{
	double rmax = 8;
	
	gsl_multiroot_fsolver *s = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_hybrid,1);

	gsl_multiroot_function F;
	F.f = master;
	F.n = 1;
	F.params = &rmax;

	gsl_vector *ev = gsl_vector_alloc(1);
	gsl_vector_set(ev, 0, -1);

	gsl_multiroot_fsolver_set(s, &F, ev);

	int status, iter = 0;

	do{
		iter++;
		status = gsl_multiroot_fsolver_iterate(s);
		if(status){
			fprintf(stderr, "Error in multiroot iterate");
			break;}
		fprintf(stderr, "Iteration %i\n", iter);

		status = gsl_multiroot_test_residual(s->f, 1e-3);
		if(status==GSL_SUCCESS)fprintf(stderr, "Function converged!\n");
	}while(status==GSL_CONTINUE && iter<100);

	double e = gsl_vector_get(s->x, 0);
	printf("%g %g\n\n\n", rmax, e);

	for(double r=0.0; r<rmax; r+=rmax/50){
		printf("%g %g %g\n", r, Fe(e,r), r*exp(-r));}
	
	gsl_multiroot_fsolver_free(s);
	gsl_vector_free(ev);
	return 0;
}





