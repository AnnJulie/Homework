#include "ODE.h"

/* Embedded midpoint-Euler method with error estimate */
void rkstep12(
        void f(int n, double x, double* yx, double* dydx),
        int n,
        double x,
        double* yx,
        double h,
        double* yh,
        double* dy)
{
	int i; 
	double k0[n], yt[n], k12[n];
	f(n, x, yx, k0);	for(i=0; i<n; i++){yt[i] = yx[i] +  k0[i]*h/2;}
	f(n, x+h/2, yt, k12);	for(i=0; i<n; i++){yh[i] = yx[i] + k12[i]*h;}

	for(i=0; i<n; i++) { dy[i] = (k0[i] - k12[i])*h/2;}
}


/* ODE driver */
int driver(
	void f(int n, double x, double* y, double* dydx),
	int n, double x0, double y0, double b, double h, double acc, double eps, int max,
	gsl_matrix* R)
{
	// allocate space for x and y values
	double* xlist = (double*)malloc(max*sizeof(double));
        double** ylist= (double**)malloc(max*sizeof(double*));
        for(int i=0; i<max; i++){ ylist[i] = (double*)malloc(n*sizeof(double)); }

	// set initial values
        xlist[0] = x0;
        ylist[0][0] = y0;
	gsl_matrix_set(R,0,0,x0);
        for(int i=0; i<n; i++) gsl_matrix_set(R,0,i+1,y0);

	// driver routine
	int i, k=0;
	double x, *y, s, err, normy, tol, a=xlist[0], yh[n], dy[n];
	while( xlist[k] < b){
		x=xlist[k], y=ylist[k]; if(x+h>b) h=b-x;
		rkstep12(f, n, x, y, h, yh, dy);
		s=0; for(i=0; i<n; i++) s+=dy[i]*dy[i]; err = sqrt(s);
		s=0; for(i=0; i<n; i++) s+=yh[i]*yh[i]; normy=sqrt(s);
		tol = (normy*eps + acc)*sqrt(h/(b-a));
		if(err<tol){ /*accept step and continue*/
			k++; if(k>max-1){
				fprintf(stderr, "Not enough allocated space for x\n");
				return -k;}
			xlist[k]=x+h;
			for(i=0; i<n; i++) {
				ylist[k][i] = yh[i];
				gsl_matrix_set(R,k,0,x+h);
                		for(i=0; i<n; i++) gsl_matrix_set(R,k,i+1,yh[i]);
			}
		}
		if(err>0) h*=pow(tol/err, 0.25)*0.95;
		else h*=2;
	}		
	// free allocated space
	for(int i=0; i<max; i++){ free(ylist[i]); }
	free(ylist);
	free(xlist);

	return k+1; /* return the number of entries in xlist/ylist */
}
	


double ODE_integral(
	void f(int n, double x, double* y, double* dydx),
        int n, double x0, double b, double h, double acc, double eps, int max,
        gsl_matrix* R)
{
	double y0 = 0;
	int last = driver(f, n, x0, y0, b, h, acc, eps, max, R);	
	return gsl_matrix_get(R, last-1, 1);
}




