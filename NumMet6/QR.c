#include "QR.h"


// Upper triangle nxn matrix R. Solves Rx = b for vector x
void backSubstitute(gsl_matrix* R, gsl_vector* x, gsl_vector* b)
{
	double hs;
	int n = R->size1-1;
	for (int i=n; i>=0; i--){
		hs = gsl_vector_get(b, i);
		for (int j=n; j>i; j--){
			hs -= gsl_matrix_get(R,i,j)*gsl_vector_get(x,j);
		} 
		gsl_vector_set(x, i, hs/gsl_matrix_get(R, i, i));
	}
}

// perform in-place modified Gram-Schmidt orthogonalization of an nxm (n>=m) matrix A.
// A is transformed to  Q and the mxm matrix R is computed.
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

// Solves the system Ax = b by QR decomposition. A and b is conserved.
void QR_solve(gsl_matrix* A, gsl_vector* x, gsl_vector* b)
{
	// allocate space
	gsl_matrix* Q = gsl_matrix_calloc(A->size1, A->size2);
	gsl_matrix* R = gsl_matrix_calloc(A->size2, A->size2);
	gsl_vector* QTb = gsl_vector_calloc(A->size2);

	// solve system
	gsl_matrix_memcpy(Q, A);
	gramschmidt(Q, R);
	gsl_blas_dgemv(CblasTrans, 1.0, Q, b, 0.0, QTb);
	backSubstitute(R, x, QTb);

	// free allocated space
	gsl_matrix_free(Q);
	gsl_matrix_free(R);
	gsl_vector_free(QTb);
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


