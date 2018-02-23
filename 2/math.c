#include "stdio.h"
#include <tgmath.h>
#include <complex.h>

int main()
{
	
	double complex z1 = csqrt(-2.0);
	double complex z2 = exp(I);
	double complex z3 = exp(I * M_PI);
	double complex z4 = pow(I,M_E);

	float xf = 0.1111111111111111111111111111f;
	double xd = 0.1111111111111111111111111111;
	long double xld = 0.11111111111111111111111L;
	
	printf("\nex1\n");
	printf("Gamma_function(5) = %f\n", tgamma(5));
	printf("Bessel function(0.5) = %f\n", j1(0.5));
	printf("sqrt(-2) = %f +  %fi\n", creal(z1), cimag(z1)); 
	printf("e^i = %f + %fi\n", creal(z2), cimag(z2));
	printf("e^(i*pi) = %f + %fi\n", creal(z3), cimag(z3));
	printf("i^e = %f + %fi\n", creal(z4), cimag(z4));

	printf("\nex2\n");
	printf("float: %.25g\ndouble: %0.25lg\nlong double: %0.25Lg\n",xf, xd, xld); 
	return 0;
}
