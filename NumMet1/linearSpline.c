#include "linearSpline.h"
#include <assert.h>


double linterp(int n, double *x, double *y, double z){
        assert(n>1 && z>x[0] && z<x[n-1]);
        int i=0, j=n-1; 
	
	while (j-i>1){ /* binary search */
		int m=(i+j)/2;if(z > x[m]) i=m;else j=m;
	}
	return y[i] + ((y[i+1]-y[i])/(x[i+1]-x[i])*(z-x[i]));
}



double linterp_integ(int n, double *x, double *y, double z){
	assert(n>1 && z>x[0] && z<x[n-1]);
        int i=0, j=n-1;

        while (j-i>1){ /* binary search */
                int m=(i+j)/2;if(z > x[m]) i=m;else j=m;
        }
	double I = 0;
	int i_max = i;
	for (i=0; i<i_max; i++){
		double h = x[i+1]-x[i];
		double a = (y[i+1]-y[i]) / h;
		I += y[i]*h + a*h*h/2;
	}
	
	double h = z-x[i_max];
	double a = (y[i_max+1]-y[i_max]) / (x[i_max+1]-x[i_max]); 
	return I + y[i]*h + a*h*h/2;
}

/*
int main(){
	int n = 25;
	double z = 2.5;
	double x[n], y[n], e1, e2;
        int i = 0;
	
 reading input table, assign to vectors, write table to out-file 
        while (scanf("%lg", &e1) != EOF && scanf("%lg", &e2) != EOF){
		x[i] = e1;
                y[i] = e2;
		printf("%g %g\n", e1, e2);
		i++;
        }

writing linear interpolation datapoint
	printf("\n\n");	
	printf("%g %g\n", z, linterp(i, x, y, z));


writing linear interpolation integral
	printf("\n\n"); Gnuplot index 2
        printf("%g %g\n", x[0], y[0]);
        printf("%g %g\n", z, linterp(i, x, y, z));
	fprintf(stderr, "\nLinear integral: I = %g\n", linterp_integ(i, x,y,z));
} */



