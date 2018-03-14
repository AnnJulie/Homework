#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>

double phi_h_phi (double x, void *params){
	double a = *(double *) params;
	double f = 0.5 * (-a*a*x*x + a + x*x)*exp(-a*x*x);
	return f;
}

double phi_phi (double x, void *params){
        double a = *(double *) params;
        double f = exp(-a*x*x);
        return f;
}

int main(){
	double a_min = .1;
	double a_max = 4;

	gsl_integration_workspace *w
		= gsl_integration_workspace_alloc(1000);

	double result_H, error_H, result_N, error_N;
	gsl_function H;
	H.function = &phi_h_phi;
	H.params = NULL;
	
	gsl_function N;
	N.function = &phi_phi;
	N.params = NULL;

	for (double a = a_min; a < a_max; a += 0.1){
		H.params = &a;
		N.params = &a;

		gsl_integration_qagi(&H, 1e-3, 1e-3, 1000, w, &result_H, &error_H);
		gsl_integration_qagi(&N, 1e-3, 1e-3, 1000, w, &result_N, &error_N);

		printf("%g %g\n", a, result_H/result_N);
	}
        gsl_integration_workspace_free(w);

        return 0;
}
