#include <stdio.h>
#include <gsl/gsl_sf_airy.h>
#include <math.h>

int main(){

	for (double x = -5; x < 5; x+=0.1){
		printf("%g %g\n", x, gsl_sf_airy_Ai(x,0));
	}

return 0;
}
