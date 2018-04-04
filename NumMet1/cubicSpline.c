#include "cubicSpline.h"
#include <stdlib.h>
#include <assert.h>

/* return index of closest point before the value z */
int binary_search_cubic(cspline *s, double z)
{
	assert (z>=s->x[0] && z<=s->x[s->n-1]);
        int i=0, j=s->n-1;
        while(j-i>1){
                int m=(i+j)/2; if(z>s->x[m]) i=m; else j=m;
        }
        return i;
}

/* allocates and builds the cubic spline */
cspline * cspline_alloc(int n, double *x, double *y)
{
	cspline* s = (cspline*)malloc(sizeof(cspline));
	s->x = (double*)malloc(n*sizeof(double));
	s->y = (double*)malloc(n*sizeof(double));
	s->b = (double*)malloc(n*sizeof(double));
	s->c = (double*)malloc((n-1)*sizeof(double));
	s->d = (double*)malloc((n-1)*sizeof(double));
	s->n = n;

	for(int i=0; i<n; i++){s->x[i]=x[i]; s->y[i]=y[i];}

	double h[n-1], p[n-1];
	for(int i=0;i<n-1;i++){h[i]=x[i+1]-x[i]; assert(h[i]>0);}
	for(int i=0;i<n-1;i++){p[i]=(y[i+1]-y[i])/h[i];}

	double D[n], Q[n-1], B[n]; //building the tridiagonal system:
	D[0]=2; D[n]=2;
	for(int i=0;i<n-2;i++){D[i+1]=2*h[i]/h[i+1]+2; }
	Q[0]=1;
	for(int i=0;i<n-2;i++){Q[i+1]=h[i]/h[i+1]; }
	B[0]=3*p[0]; B[n-1]=3*p[n-2];
	for(int i=0;i<n-2;i++){B[i+1]=3*(p[i]+p[i+1]*h[i]/h[i+1]); }
	for(int i=1;i<n;i++){ D[i]-=Q[i-1]/D[i-1]; B[i]-=B[i-1]/D[i-1]; }
	s->b[n-1]=B[n-1]/D[n-1];
	for(int i=n-2;i>0;i--){ s->b[i]=(B[i]-Q[i]*s->b[i+1])/D[i]; }
	for(int i=0;i<n-1;i++){
		s->c[i] = (-2*s->b[i] - s->b[i+1] + 3*p[i]) / h[i];
		s->d[i] = (s->b[i] + s->b[i+1] - 2*p[i]) / h[i] / h[i];
	}
	return s;
}



/* evaluates the prebuilt spline at point z */
double cspline_evaluate(cspline *s, double z)
{
	int i = binary_search_cubic(s, z);
        double h=z-s->x[i];

        return s->y[i] + h*(s->b[i] + h*(s->c[i] + h*s->d[i]));
}


/* evaluates the derivative of the prebuilt spline at point z */
double cspline_derivative(cspline *s, double z)
{
	int i = binary_search_cubic(s, z);
        double h = z-s->x[i];

        return s->b[i] + 2*s->c[i]*h + 3*s->d[i]*h*h;
}

/* evaluates the integral of the prebuilt spline from x[0] to z */
double cspline_integral(cspline *s, double z)
{
	double I = 0;
        int i_max = binary_search_cubic(s, z);

        for (int i=0;i<i_max;i++){
                double h = s->x[i+1] - s->x[i];
                I += s->y[i]*h + s->b[i]*h*h/2 + s->c[i]*h*h*h/3 + s->d[i]*h*h*h*h/4;
        }

        double h = z-s->x[i_max];
        return I + s->y[i_max]*h + s->b[i_max]*h*h/2 + s->c[i_max]*h*h*h/3 + s->d[i_max]*h*h*h*h/4;
}

/* free memory allocated in qspline_alloc */
void cspline_free(cspline *s)
{
	free(s->x); free(s->y); free(s->b); free(s->c); free(s->d); free(s);
}
