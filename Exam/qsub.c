#include "qsub.h"

qsub_spline* qsub_spline_alloc(int n, double *x, double *y, double *dy)
{
	assert(n>2);
	double h[n-1], p[n-1];
	for(int i=0; i<n-1; i++){h[i]=x[i+1]-x[i]; assert(h[i]>0);}
	for(int i=0; i<n-1; i++){p[i]=(y[i+1]-y[i])/h[i];}

	qsub_spline* s = (qsub_spline*)malloc(sizeof(qsub_spline));
	s->x = (double*)malloc(n*sizeof(double));
	s->y = (double*)malloc(n*sizeof(double));
	s->dy = (double*)malloc(n*sizeof(double));
	s->c = (double*)malloc((n-1)*sizeof(double));
	s->d = (double*)malloc((n-1)*sizeof(double));
	s->n = n;
	for(int i=0; i<n; i++){s->x[i]=x[i]; s->y[i]=y[i]; s->dy[i]=dy[i];}

	for(int i=0; i<n-1; i++){
		s->c[i] = (3*p[i] - 2*s->dy[i] - s->dy[i+1])/h[i];
		s->d[i] = (s->dy[i+1] + s->dy[i] - 2*p[i])/(h[i]*h[i]);
	}

	return s;
}


double qsub_spline_eval(qsub_spline* s, double z)
{
	assert(z>s->x[0] && z<s->x[s->n-1]);
	int i=0, j=s->n-1;
	while(j-i > 1){int m=(i+j)/2; if(z>s->x[m]) i=m; else j=m;}	
	double h = z - s->x[i];
	return s->y[i] + h*(s->dy[i] + h*(s->c[i] + h*s->d[i]));
}


void qsub_spline_free(qsub_spline* s)
{
	free(s->x);
	free(s->y);
	free(s->dy);
	free(s->c);
	free(s->d);
	free(s);
}




