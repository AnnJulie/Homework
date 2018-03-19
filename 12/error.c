#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>


int ode(double x, const double u[], double dudx[], void *params)
{
	(void)(x); /* avoid unused parameter warning */

	dudx[0] = 2/sqrt(M_PI) * exp(-x*x);
	return GSL_SUCCESS;
}



int main(int argc, char *argv[])
{
	/* reading parameters from Makefile-input */
	double a = atof(argv[1]);
	double b = atof(argv[2]);
	double dx = atof(argv[3]);
	fprintf(stderr, "\nParameters:\na = %g\nb = %g\ndx = %g\n\n", a, b, dx);

	/* setting up environment */
	gsl_odeiv2_system sys = {ode, NULL, 1, NULL};
	gsl_odeiv2_driver *driver = 
		gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0);

	/* solving the system */
	double t=0;
	double y[1] = {0.0};

	/* data for x>0 */
	for (double x=0; x<b; x+=dx)
	{
		gsl_odeiv2_driver_apply(driver, &t, x, y);
		printf("%g %g\n", x, y[0]);
	}
	printf("\n");

	/* data for x<0 */
	t=0;
	y[0]=0;
	for (double x=0; x<-a; x+=dx)
        {
                gsl_odeiv2_driver_apply(driver, &t, x, y);
                printf("%g %g\n", -x, -y[0]);
        }

	gsl_odeiv2_driver_free(driver);
	return 0;
}





	
