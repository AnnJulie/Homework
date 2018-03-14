#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>

int diff_ode(double x, const double y[], double dydx[], void *params)
{
	(void)(x); /* avoid unused parameter warning */
	
	dydx[0] = y[0]*(1.0-y[0]);
	return GSL_SUCCESS;
}


double my_diff(double x){
	gsl_odeiv2_system sys = {diff_ode, NULL, 1, NULL}; 


	gsl_odeiv2_driver *driver = 
		gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0);


	double t=0;
	double y[1] = {0.5};
	gsl_odeiv2_driver_apply(driver, &t, x, y);

	gsl_odeiv2_driver_free(driver);
	return y[0];
}


int main(){
	for (int x = 0; x < 300; x++){
		printf("%g, %g \n", x/100.0, my_diff(x/100.0));
	}

	return 0;
}
