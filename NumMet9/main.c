#include "MC.h"
#include "f.h"
#include <math.h>
#include <stdio.h>
#define PI 3.14159265358979323846
int main()
{
	printf("Integrating f(x,y) = x*y from x=[0,2], y=[0,2]:\n");
	int N1;
	int dim1 = 2;
	double a1[dim1];
	double b1[dim1];
	double* result1 = (double*)malloc(dim1*sizeof(double));
	double* error1 = (double*)malloc(dim1*sizeof(double));
	a1[0] = 0;
	a1[1] = 0;
	b1[0] = 2;
	b1[1] = 2;
        printf("Using N = 10\n");
        N1 = 10;
        plainmc(dim1, a1, b1, f1, N1, result1, error1);
	printf("Using N = 100\n");
	N1 = 100;
	plainmc(dim1, a1, b1, f1, N1, result1, error1);
        printf("Using N = 1000\n");
        N1 = 1000;
        plainmc(dim1, a1, b1, f1, N1, result1, error1);
        printf("Using N = 10000\n");
        N1 = 10000;
        plainmc(dim1, a1, b1, f1, N1, result1, error1);
	printf("Increasing N by a factor of 10 results in a decreace in the error of approximately sqrt(10) ≈ 3.16\n");

	printf("\n\nIntegrating f(x,y,z) = (1/π^3)*[1-cos(x)cos(y)cos(z)]^-1 from x=[0,pi], y=[0,pi], z=[0,pi]:\n");
        int dim2 = 3;
        int N2 = 1000000;
        double a2[dim2];
        double b2[dim2];
        double* result2 = (double*)malloc(dim2*sizeof(double));
        double* error2 = (double*)malloc(dim2*sizeof(double));
        a2[0] = 0;
        a2[1] = 0;
	a2[2] = 0;
        b2[0] = PI;
        b2[1] = PI;
	b2[2] = PI;
        plainmc(dim2, a2, b2, f2, N2, result2, error2);
	
	return 0;
}

