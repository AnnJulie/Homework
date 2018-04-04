#include "gramSchmidt.h"
#include "matrix.h"
#include "vector.h"
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_blas.h>

int main()
{
	// generate random 3x2 matrix A and vector b
	gsl_matrix* A = gsl_matrix_calloc(3,2);
	gsl_vector* b = gsl_vector_calloc(3);
//	matrix_generate_random(A);
	gsl_matrix_set(A, 0, 0, 3);
	gsl_matrix_set(A, 1, 0, 4);
	gsl_matrix_set(A, 2, 0, 0);
	gsl_matrix_set(A, 0, 1, -6);
	gsl_matrix_set(A, 1, 1, -8);
	gsl_matrix_set(A, 2, 1, 1);
        gsl_vector_set(b,0,3);
        gsl_vector_set(b,1,4);
        gsl_vector_set(b,2,2);
        printf("Generate matrix A and vector b for equation Ax = b\n");
        vector_print(b, "b");
	matrix_print(A, "A");

	// create zero matrix R
	gsl_matrix* R = gsl_matrix_calloc(2,2);

	// perform Gram-Schmidt factorization
	qr_gs_decomp(A, R);
	printf("Perform QR decomposition\n");	
	matrix_print(A, "Q");
	matrix_print(R, "R");


	// testing solution: QR = A; (Q^T)Q = 1, where Q is stored in the matrix A after
	// Gram-Schmidt factorization is performed
	gsl_matrix* QR = gsl_matrix_calloc(3,2);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A, R, 0.0, QR); 
	printf("Check QR = A:\n");
	matrix_print(QR, "QR");

	printf("Check QTQ = 1:\n");
	gsl_matrix* QTQ = gsl_matrix_calloc(2,2);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, A, A, 0.0, QTQ);
	matrix_print(QTQ, "(Q^T)Q");


	gsl_vector* x = gsl_vector_calloc(2);
	qr_gs_solve(A, R, b, x);
	printf("Solution!\n");
	vector_print(x, "x");
	
	// Check solution: QRx = b
	gsl_vector* QRx = gsl_vector_calloc(3);
	gsl_blas_dgemv(CblasNoTrans, 1.0, QR, x, 0.0, QRx);
	printf("Checking solution:\n");
	vector_print(QRx, "Ax (should be equal to b)");


	// Finding inverse of A:
	printf("Compute B = A-1\n");
	gsl_matrix* B = gsl_matrix_calloc(2,3);
	qr_gs_inverse(A, R, B);
	matrix_print(B, "B");

	printf("Check BA = 1\n");
	gsl_matrix* AB = gsl_matrix_calloc(2,2);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, B, QR, 0.0, AB);
	matrix_print(AB, "BA");

	// free allocated space
	gsl_matrix_free(A);
	gsl_matrix_free(R);
	gsl_matrix_free(QR);
	gsl_matrix_free(QTQ);
	gsl_matrix_free(B);
	gsl_matrix_free(AB);
	gsl_vector_free(b);
	gsl_vector_free(x);
	gsl_vector_free(QRx);
}
