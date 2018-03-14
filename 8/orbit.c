#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <unistd.h>


int orbit_ode(double phi, const double u[], double dudphi[], void *params)
{
	(void)(phi); /*avoid unused parameters*/

	double epsilon = *(double *) params;	

	dudphi[0] = u[1];
	dudphi[1] = 1-u[0] + epsilon * u[0] * u[0];
	return GSL_SUCCESS;
}



int main(int argc, char *argv[])
{
	
	double epsilon = 0, dudphi=0;
	/* reading parameters from short command-line optons */
	while(1){
		int opt = getopt(argc, argv, "e:p:");
		if(opt == -1) break;
		switch(opt){
			case 'e': epsilon=atof(optarg);break;
			case 'p': dudphi=atof(optarg); break;
			default:
				printf("failiure setting parameters");
				exit(EXIT_FAILURE);
		}
	};
	

	gsl_odeiv2_system sys ={orbit_ode, NULL, 2, (void *) &epsilon};

	gsl_odeiv2_driver *driver =
		gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, 1e-3, 1e-6, 1e-6);

	double t = 0;
	double u[2] = {1, dudphi};
	double phi_max = 95 * M_PI, delta_phi = 0.05;	

	for (double phi = 0; phi < phi_max; phi += delta_phi){

		int status = gsl_odeiv2_driver_apply(driver, &t, phi, u);
	
		printf("%g %g\n", phi, u[0]);
	}

	gsl_odeiv2_driver_free (driver);
	return 0;
}	




