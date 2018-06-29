#include "qsub.h"
#include <stdio.h>
#include <math.h>

int main()
{
	/* reading input-values to n, x_min and x_max */
	int n;
	double x_min, x_max, step_size;
        scanf("%d",&n);
        scanf("%lg",&x_min);
        scanf("%lg",&x_max);
	scanf("%lg",&step_size);

	/* reading input values to x[n], y[n] and dy[n]. Write data-points (x,y) to out-file */
        double x[n], y[n], dy[n], e1, e2, e3;
        int i = 0;
        while (scanf("%lg",&e1)!=EOF && scanf("%lg",&e2)!=EOF && scanf("%lg",&e3)!=EOF &&  i<n){
                x[i] = e1;
                y[i] = e2;
		dy[i]= e3;
                printf("%g %g\n", e1, e2);
                i++;
        }

	/* Initializing the qubic sub-spline */
	qsub_spline* s = qsub_spline_alloc(n, x, y, dy);	

        /* printing evaluation from the qubic sub-spline */
        printf("\n\n"); /* Gnuplot index 1 */
	double z = x_min + step_size;
	while(z<x_max){
	        printf("%g %g\n", z, qsub_spline_eval(s, z));
		z += step_size;
	}
	
}
