#include <stdio.h>
#include <math.h>

int main(){

	for (int x = 0; x < 30; x++){
		printf("%g %g, \n", x/10.0, exp(x/10.0)/(1+exp(x/10.0)) );
	}

	return 0;
}
