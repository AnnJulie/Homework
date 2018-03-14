#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>

double f (double x, void *params){

	double f = log(x) / sqrt(x);
	return f;
}

int main (){

	gsl_integration_workspace *w
		= gsl_integration_workspace_alloc(1000);

	double result, error;
	gsl_function F;
	F.function = &f;
	F.params = NULL;

	gsl_integration_qags(&F, 0, 1, 1e-3, 1e-3, 1000, w, &result, &error);

	printf("result = % .18f\n", result);
	printf("estimated error = % .18f\n", error);

	gsl_integration_workspace_free(w);

	return 0;
}
