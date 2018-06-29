#include "f.h"
#include "MC.h"
#include <omp.h>
#include <math.h>
#include <time.h>
#define PI 3.14159265358979323846

int main()
{
	int N = 10000000;
	int dim = 1;
	printf("TASK:\n");
	printf("Integrate a function defined differently overt three intervals\n");
	printf("          x*cos(x)       ,    0 < x < PI\n");
	printf("  f(x) =  sin(x)*sin(x)  ,   PI < x < 2*PI\n");
	printf("          cos(x)*cos(x)  , 2*PI < x < 3*PI\n\n");
        double a1[dim], b1[dim], a2[dim], b2[dim], a3[dim], b3[dim];
        double* result1 = (double*)malloc(dim*sizeof(double));
        double* error1 = (double*)malloc(dim*sizeof(double));
	double* result2 = (double*)malloc(dim*sizeof(double));
        double* error2 = (double*)malloc(dim*sizeof(double));
        double* result3 = (double*)malloc(dim*sizeof(double));
        double* error3 = (double*)malloc(dim*sizeof(double));

        a1[0] = 0;
        b1[0] = PI;
        a2[0] = PI;
        b2[0] = 2*PI;
	a3[0] = 2*PI;
	b3[0] = 3*PI;

	printf("first integrate the intervals after each other:\n");
	clock_t begin = clock();
	plainmc(dim, a1, b1, f1, N, result1, error1);
	plainmc(dim, a2, b2, f2, N, result2, error2);
	plainmc(dim, a3, b3, f3, N, result3, error3);
	printf("result = %g, error = %g\n", 
                *result1 + *result2 + *result3,
                *error1 + *error2 + *error3);
	clock_t end = clock();
	double time = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("time = %g sec\n", time);



	printf("\nNow integrating the intervals in  parallel sections:\n");
	begin = clock();
	#pragma omp parallel sections
	{
		#pragma omp section
		{
        		plainmc(dim, a1, b1, f1, N, result1, error1);

		}
		#pragma omp section
		{
                        plainmc(dim, a2, b2, f2, N, result2, error2);
		}

		#pragma omp section
                {
                        plainmc(dim, a3, b3, f3, N, result3, error3);
                }
	}
	printf("result = %g, error = %g\n", 
		*result1 + *result2 + *result3, 
		*error1 + *error2 + *error3);
	end = clock();
        time = (double)(end - begin) / CLOCKS_PER_SEC;
        printf("time = %g sec\n\n", time);


	printf("It took longer time to split the processor in three parts than the gain of integrating simultaneously. This thus serves as a bad example of how to save time by multiprocessing\n");
}



