#include "stdio.h"
#include "math.h"

int equal(double a, double b, double tau, double epsilon)
{
	printf("\na = %f\n", a);
	printf("b = %f\n", b);
	printf("tau = %f\n", tau);
	printf("epsilon = %f\n\n", epsilon);


	if (fabs(a-b) < tau){
		printf("|a-b| < tau:\nreturning 1\n");
		return 1;
	}else if(fabs(a-b)/(fabs(a)+fabs(b)) < epsilon){
		printf("|a-b|/(|a|+|b|) < epsilon:\nreturning 1\n");
	}else{
		printf("Numbers are not equal:\nreturning 0\n");
		return 0;
	}

	// just for safety
	return 0;
}
