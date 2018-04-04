#include "quadraticSpline.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
 

int binary_search(qspline *s, double z)
{
	assert (z>=s->x[0] && z<=s->x[s->n-1]);
        int i=0, j=s->n-1;
        while(j-i>1){
                int m=(i+j)/2; if(z>s->x[m]) i=m; else j=m;
        }
	return i;
}

void vector_print(int n, double *v)
{
	printf("[");
	for (int i=0; i<n-1; i++){
		printf("%g,",v[i]);
	}
	printf("%g]\n", v[n]);
}



/* allocates and builds the quadratic spline */
qspline* qspline_alloc(int n, double *x, double *y, int print)
{
	qspline *s = (qspline*)malloc(sizeof(qspline));
	s->b = (double*)malloc((n-1)*sizeof(double));
	s->c = (double*)malloc((n-1)*sizeof(double));
	s->x = (double*)malloc(n*sizeof(double));
	s->y = (double*)malloc(n*sizeof(double));	
	s->n = n;
	
	for (int i=0;i<n;i++){
		s->x[i] = x[i];
		s->y[i] = y[i];
	}

	int i; double p[n-1], h[n-1];
	for (i=0;i<n-1;i++){
		h[i] = x[i+1]-x[i];
		p[i] = (y[i+1]-y[i])/h[i];
	}
	
	s->c[0] = 0;
	for (i=0;i<n-2;i++){ //recursion up
		s->c[i+1] = (p[i+1]-p[i] - s->c[i]*h[i]) / h[i+1];
	}

	s->c[n-2]/=2; // recursion down
	for (i=n-3;i>=0;i--){
		s->c[i] = (p[i+1]-p[i] - s->c[i+1]*h[i+1]) / h[i];
	}

	for (i=0;i<n-1;i++){
		s->b[i] = p[i]- s->c[i]*h[i];
	}
	
	if(print){
		printf("Printing vector b\n");
		vector_print(n-1, s->b);
		printf("Printing vector c\n");
		vector_print(n-1, s->c);
	}
	return s;
}


/* evaluates the prebuilt spline at point z */
double qspline_evaluate(qspline *s, double z)
{
	int i = binary_search(s, z);
	double h=z-s->x[i];

	return s->y[i] + h*(s->b[i] + h*s->c[i]); // interpolating polynomial
}


/* evaluates the derivative of the prebuilt spline at point z */
double qspline_derivative(qspline *s, double z)
{
	int i = binary_search(s, z);
	double h = z-s->x[i];

	return s->b[i] + 2*s->c[i]*h;
}


/* evaluates the integral of the prebuilt spline from x[0] to z */
double qspline_integral(qspline *s, double z)
{
	double I = 0;
	int i_max = binary_search(s, z);

	for (int i=0;i<i_max;i++){
		double h = s->x[i+1] - s->x[i];
		I += s->y[i]*h + s->b[i]*h*h/2 + s->c[i]*h*h*h/3; 
	}
	
	double h = z-s->x[i_max];
	return I + s->y[i_max]*h + s->b[i_max]*h*h/2 + s->c[i_max]*h*h*h/3;
}

/* free memory allocated in qspline_alloc */
void qspline_free(qspline *s){
	free(s->x); free(s->y); free(s->b); free(s->c); free(s);
}
