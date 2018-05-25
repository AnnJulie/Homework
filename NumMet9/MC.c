#include "MC.h"
#define RND ((double)rand()/RAND_MAX)

void randomx(int dim, double* a, double* b, double* x)
{
	for (int i=0; i<dim; i++) x[i] = a[i] + RND*(b[i]-a[i]);
}


void plainmc(int dim, double* a, double* b, double f(double* x), int N, double* result, double* error)
{
	double V=1;
	for(int i=0; i<dim; i++) V *= b[i]-a[i];

	double sum=0, sum2=0, fx, x[dim];
	for(int i=0; i<N; i++){
		randomx(dim, a, b, x);
		fx = f(x);
		sum += fx;
		sum2 += fx*fx; }
	
	double avr = sum/N, var = sum2/N - avr*avr;
	*result = avr*V;
	*error = sqrt(var/N)*V;
	printf("result = %g, error = %g\n", *result, *error);
	return;
}
