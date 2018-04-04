#include "gramSchmidt.h"


void qr_gs_inverse(const gsl_matrix* Q, const gsl_matrix* R, gsl_matrix* B)
{
	// finding inverse of R
        int s;
        gsl_permutation* p = gsl_permutation_alloc(2);
        gsl_matrix* invR = gsl_matrix_calloc(R->size1, R->size2);
        gsl_linalg_LU_decomp(R, p, &s);
        gsl_linalg_LU_invert(R, p, invR);

	// inverse of A=QR is invR(Q^T)
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, invR, Q, 0.0, B);

	gsl_permutation_free(p);
	gsl_matrix_free(invR);
}




// solving equation QRx = b for x
void qr_gs_solve(const gsl_matrix* Q, gsl_matrix* R,const gsl_vector* b, gsl_vector* x){
	// performimg (Q^T)b, storing solution in b_temp
	gsl_vector* QTb = gsl_vector_calloc(2);
	gsl_blas_dgemv(CblasTrans, 1.0, Q, b, 0.0, QTb);
	
	// finding inverse of R
	int s;
	gsl_permutation* p = gsl_permutation_alloc(2);
	gsl_matrix* invR = gsl_matrix_calloc(R->size1, R->size2);
	gsl_linalg_LU_decomp(R, p, &s);
	gsl_linalg_LU_invert(R, p, invR);

	// performing invR(Q^T)b_temp, storing solution in x
	gsl_blas_dgemv(CblasNoTrans, 1.0, invR, QTb, 0.0, x);
	
	// free allocated space
	gsl_permutation_free(p);
	gsl_vector_free(QTb);
	gsl_matrix_free(invR);
}


// perform in-place modified Gram-Schmidt orthogonalization of an nxm (n>=m) matrix A.
// A turns into Q and the square mxm matrix R is computed.
void qr_gs_decomp(gsl_matrix *A, gsl_matrix* R)
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

