#include "integrator.h"


double adapt24(double f(double, gsl_vector* calls), double a, double b, double acc, double eps, double f2, double f3, int nrec, gsl_vector* calls)
{
	assert(nrec<100000000);
	double f1 = f(a+(b-a)/6, calls);
	double f4 = f(a+5*(b-a)/6, calls);
	double Q = (2*f1+f2+f3+2*f4)/6*(b-a);
	double q = (f1+f2+f3+f4)/4*(b-a);
	double tolerance = acc+eps*fabs(Q);
	double error = fabs(Q-q);
	
	if(error < tolerance) return Q;
	
	double Q1 = adapt24(f,a,(a+b)/2,acc/sqrt(2.),eps,f1,f2,nrec+1,calls);
	double Q2 = adapt24(f,(a+b)/2,b,acc/sqrt(2.),eps,f3,f4,nrec+1,calls);
	return Q1+Q2;
}


double adapt(double f(double, gsl_vector* calls), double a, double b, double acc, double eps, gsl_vector* calls)
{
	if(isinf(a) != 0){a=-2100;}
	if(isinf(b) != 0){b=2100;}

	gsl_vector_set(calls, 0, 1);
	double f2 = f(a+2*(b-a)/6, calls);
	double f3 = f(a+4*(b-a)/6, calls);
	int nrec = 0;
	return adapt24(f,a,b,acc,eps,f2,f3,nrec,calls);
}









