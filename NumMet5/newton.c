#include "newton.h"


void newton(void f(gsl_vector* x, gsl_vector* fx), gsl_vector* xstart)
{
	double eps = 1e-3; // error
	double dx = 1e-6;  // step
	double l = 2;	   // lambda
	size_t n = xstart->size;
	int count = 0;	   // safety

	// allocate space
	gsl_matrix* J = gsl_matrix_calloc(n,n);
	gsl_vector* fx = gsl_vector_calloc(n);
	gsl_vector* df = gsl_vector_calloc(n);
	gsl_vector* fy = gsl_vector_calloc(n);
	gsl_vector* Dx = gsl_vector_calloc(n);
	gsl_vector* y = gsl_vector_calloc(n);
	gsl_vector* b = gsl_vector_calloc(n);
	gsl_vector* x = gsl_vector_calloc(n);

	// initial values
	gsl_vector_memcpy(x, xstart);
	f(x, fx);
	
	// iterative method
	while(count<100){
	count++;
	for (int j=0; j<n; j++){
		gsl_vector_set(x, j, gsl_vector_get(x,j)+dx);
		f(x, df); //  df is now equal to f(x+dx)
		gsl_vector_sub(df, fx); // df is now equal to f(x+dx)-f(x)
		for (int i=0; i<n; i++){
			gsl_matrix_set(J, i, j, gsl_vector_get(df, i)/dx);
		}
		gsl_vector_set(x, j, gsl_vector_get(x,j)-dx);
	}

	// Dx = solution to Jx=-fx
	gsl_vector_memcpy(b, fx);
	gsl_vector_scale(b, -1);
	QR_solve(J, Dx, b);	

	// backtracking to find lambda
	l = 2;
	while (1){
		l/=2;
		gsl_vector_memcpy(y, x);  // y = x
		gsl_blas_daxpy(l, Dx, y); // y = l*Dx + y
		f(y, fy);
		if(gsl_blas_dnrm2(fy) < (1-l/2)*gsl_blas_dnrm2(fx) || l < 0.02){break;}
	}
	gsl_vector_memcpy(x, y);
	gsl_vector_memcpy(fx, fy);

	// convergence criteria
	if(gsl_blas_dnrm2(Dx) < dx || gsl_blas_dnrm2(fx) < eps){
		printf("Converged after %d iterations.\n", count);
		vector_print(x, "x");
		break;
	}
	// safety
	if(count==99){printf("count = %d\n", count);}
	}


	// free allocated space
	gsl_matrix_free(J);
	gsl_vector_free(fx);
	gsl_vector_free(df);
	gsl_vector_free(fy);
	gsl_vector_free(Dx);
	gsl_vector_free(y);
	gsl_vector_free(b);
	gsl_vector_free(x);
}


void newton_with_jacobian(void f(gsl_vector* x, gsl_vector* fx, gsl_matrix* J), gsl_vector* xstart)
{
	double eps = 1e-3; // error
	double l = 2;	   // lambda
	int count = 0;	   // safety
	size_t n = xstart->size;

	// allocation
	gsl_matrix* J = gsl_matrix_calloc(n,n); 
	gsl_vector* fx = gsl_vector_calloc(n);  
	gsl_vector* fy = gsl_vector_calloc(n);  
	gsl_vector* Dx = gsl_vector_calloc(n); 
	gsl_vector* y = gsl_vector_calloc(n);   
	gsl_vector* b = gsl_vector_calloc(n);   
	gsl_vector* x = gsl_vector_calloc(n);

  	// initial values
	gsl_vector_memcpy(x, xstart);
	f(x, fx, J);

	// iterative method
	while(count<100){
	count++;

	// Dx = solution to Jx=-fx
	gsl_vector_memcpy(b, fx);
	gsl_vector_scale(b, -1);
	QR_solve(J, Dx, b);	

	// backtracking to find lambda
	l = 2;
	while (1){
		l/=2;
		gsl_vector_memcpy(y, x);  // y = x
		gsl_blas_daxpy(l, Dx, y); // y = s*Dx + y
		f(y, fy, J);
		if(gsl_blas_dnrm2(fy) < (1-l/2)*gsl_blas_dnrm2(fx) || l < 0.02){break;}
	}
	gsl_vector_memcpy(x, y);
	gsl_vector_memcpy(fx, fy);

	// convergence criteria
	if(gsl_blas_dnrm2(fx) < eps){
		printf("Converged after %d iterations.\n", count);
		vector_print(x, "x");
		break;
	}
	// safety
	if(count==99){printf("count = %d\n", count);}
	}


	// free allocated space
	gsl_matrix_free(J);
	gsl_vector_free(fx);
	gsl_vector_free(fy);
	gsl_vector_free(Dx);
	gsl_vector_free(y);
	gsl_vector_free(b);
	gsl_vector_free(x);
}

