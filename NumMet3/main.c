#include "jacobiDiagonalization.h"
#include "matrix.h"
#include "vector.h"
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <stdio.h>



int main()
{
	testCyclicDiagonalization();
	size_t n = 12;


	printf("Testing diagonalization times:\n\n");
	printf("CYCLIC\n");
	printf("%zux%zu matrix: \ntime = %g sec\n\n", 1*n, 1*n, cyclicDiagTime(1*n));
	printf("%zux%zu matrix: \ntime = %g sec\n\n", 2*n, 2*n, cyclicDiagTime(2*n));
	printf("%zux%zu matrix: \ntime = %g sec\n\n", 3*n, 3*n, cyclicDiagTime(3*n));
	printf("n doubles and time gets approx 8 times as large. -> ~O(n)^3\n\n");

	printf("VALUE-BY-VALUE\n");
	printf("%zux%zu matrix: \ntime = %g sec\n\n", 1*n, 1*n, valueDiagTime(1*n));
	printf("%zux%zu matrix: \ntime = %g sec\n\n", 2*n, 2*n, valueDiagTime(2*n));
	printf("%zux%zu matrix: \ntime = %g sec\n\n", 3*n, 3*n, valueDiagTime(3*n));
	printf("Value-by-value method takes longer time.\n");

}
