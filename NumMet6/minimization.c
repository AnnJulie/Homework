#include "minimization.h"



/* f: objective function, df: gradient, H: Hessian matrix H*/
void newton(
	double f(gsl_vector* x, gsl_vector* df, gsl_matrix* H),
	gsl_vector* x)
{
	double fx, innerp; // storage
	double eps = 1e-3; // error
	double l = -0.2;   // lambda
	double a = 0.0001; // alpha
	size_t n = x->size;// size
	int count = 0;	   // safety

	// allocate memory
	gsl_matrix* H = gsl_matrix_calloc(n,n);
	gsl_vector* df = gsl_vector_calloc(n);
	gsl_vector* dx = gsl_vector_calloc(n);
	gsl_vector* x_temp = gsl_vector_calloc(n);

	fx = f(x, df, H); // Setting  H(x), df(x) and f(x)
	
	while(count < 30){
		count++;
		QR_solve(H, dx, df); // Solving H(x) ∆x = df(x) for ∆x

		// backtracking to find lambda
		l = -2;
		gsl_blas_ddot(dx, df, &innerp); // saving ∆xT*df(x)
		while (1){
			l/=2;
			gsl_vector_memcpy(x_temp, x);  // "remember" the value of x
			gsl_blas_daxpy(l, dx, x_temp); // x_temp = l*∆x + x_temp
			if(f(x_temp, df, H) < fx + a*l*innerp || l < 0.02){break;}
			}

		gsl_blas_daxpy(l, dx, x); // x = l*∆x + x

		// convergence criteria
		if(gsl_blas_dnrm2(df) < eps){
			printf("Converged after %d iterations.\n", count);
			vector_print(x, "x");
			break;
		}
				
		// safety
		if(count==99){printf("count = %d\n", count);}
	}

	// free allocated space
	gsl_matrix_free(H);
	gsl_vector_free(df);
	gsl_vector_free(dx);
	gsl_vector_free(x_temp);
}




