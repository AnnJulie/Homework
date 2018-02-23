#include "stdio.h"
#include <math.h>

int main()
{
	int safety = 0;
	double x;
	while( scanf("%lg", &x) != EOF ) {
		printf("%lg \t %lg\n", x, cos(x));
		safety++;
		if (safety > 10){return 0;}
	}
	
	return 0;
}
