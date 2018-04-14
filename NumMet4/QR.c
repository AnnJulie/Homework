#include "QR.h"

// allocates and build a structure to hold the data to be fitted
fit_data* fit_data_alloc(int n, double* x, double* y, double* dy){
	fit_data* data = (fit_data*)malloc(sizeof(fit_data));
	data->x = (double*)malloc(n*sizeof(double));
	data->y = (double*)malloc(n*sizeof(double));
	data->dy= (double*)malloc(n*sizeof(double));
	data->n = n;
	for(int i=0; i<n; i++){
		data->x[i] = x[i];
		data->y[i] = y[i];
		data->dy[i]= dy[i];
	}
	return data;
} 

void fit_data_free(fit_data* data)
{
	free(data->x); free(data->y); free(data->dy); free(data);
}

double funs(int i, double x)
{
	switch(i){
		case 0: return 1; break;
		case 1: return x; break;
		case 2: return x*x; break;
		default: {fprintf(stderr, "funs: wrong i: %d", i); return NAN;}
	}
}


void lsfit(int fs, fit_data* data, gsl_vector* c, gsl_matrix* S)
{
	int n = data->n;
	int m = fs;
	gsl_matrix* A = gsl_matrix_calloc(n,m);
	gsl_vector* b = gsl_vector_calloc(n);

	for (int i=0; i<n; i++){
		gsl_vector_set(b, i, data->y[i] / data->dy[i]);
		for (int k=0; k<m; k++){
			double d = data->x[i] / data->dy[i];
			gsl_matrix_set(A, i, k, funs(k, d));
		}
	}

	gsl_matrix* R = gsl_matrix_calloc(m,m);
	gramschmidt(A, R); // A is transformed to Q

	gsl_vector* QTb = gsl_vector_calloc(m);	
	gsl_blas_dgemv(CblasTrans, 1.0, A, b, 0.0, QTb);
	backSubstitute(R, c, QTb);

	gsl_matrix* invS = gsl_matrix_calloc(m,m);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, R, R, 0.0, invS);
	matrix_inverse(invS, S);

	// free allocated space
	gsl_matrix_free(A);
	gsl_matrix_free(R);
	gsl_matrix_free(invS);
	gsl_vector_free(b);
	gsl_vector_free(QTb);
}

// Upper triangle nxn matrix A. Solves Ax = b for vector x
void backSubstitute(gsl_matrix* A, gsl_vector* x, gsl_vector* b)
{
	double hs;
	int n = A->size1-1;
	for (int i=n; i>=0; i--){
		hs = gsl_vector_get(b, i);
		for (int j=n; j>i; j--){
			hs -= gsl_matrix_get(A,i,j)*gsl_vector_get(x,j);
		} 
		gsl_vector_set(x, i, hs/gsl_matrix_get(A, i, i));
	}
}

// perform in-place modified Gram-Schmidt orthogonalization of an nxm (n>=m) matrix A.
// A turns into Q and the square mxm matrix R is computed.
void gramschmidt(gsl_matrix *A, gsl_matrix* R)
{
	size_t n = A->size1;
	size_t m = A->size2;

	gsl_vector* a = gsl_vector_calloc(n);
	gsl_vector* q = gsl_vector_calloc(n);
	
	for (int i=0; i<m; i++){
		vector_get_from_matrix(a, A, i);
		gsl_matrix_set(R, i, i, sqrt(vector_inner_product(a, a)));
		vector_store_normalized(a, q);
		vector_set_to_matrix(q, A, i);
		for (int j=i+1; j<m; j++){
			vector_get_from_matrix(a, A, j);
			gsl_matrix_set(R, i, j, vector_inner_product(a, q));
			vector_update_orthogonalized(a, q);
			vector_set_to_matrix(a, A, j);
		}
	}
	gsl_vector_free(a);
	gsl_vector_free(q);
}


// generates a 3x2 matrix and vector to use in QR decomposition
void QR_generate_example(gsl_matrix* A, gsl_vector* b)
{
        gsl_matrix_set(A, 0, 0, 3);
        gsl_matrix_set(A, 1, 0, 4);
        gsl_matrix_set(A, 2, 0, 0);
        gsl_matrix_set(A, 0, 1, -6);
        gsl_matrix_set(A, 1, 1, -8);
        gsl_matrix_set(A, 2, 1, 1);

        gsl_vector_set(b,0,3);
        gsl_vector_set(b,1,4);
        gsl_vector_set(b,2,2);
}


