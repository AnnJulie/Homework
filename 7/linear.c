#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <stdio.h>
#include <stdlib.h>
#define RND (double)rand()/RAND_MAX

void matrix_print(const char* s, const gsl_matrix* m){
	printf("%s\n",s);
	for (int j=0; j<(*m).size2; j++){
		for (int i=0; i<(*m).size1; i++){
			printf("%8.3g ",gsl_matrix_get(m,i,j));
		}
		printf("\n");
	}
	printf("\n");
}


void vector_print(const char* s, const gsl_vector* v){
	printf("%s\n",s);
	for (int i=0; i<(*v).size; i++){
		printf("%8.3g ",gsl_vector_get(v,i));
	}
	printf("\n\n");
}



int main(){
	printf("Looking at the system Mx = b:\n\n");
	
	int n=3;
	int i=0;
	int j=0;
	double element;
	gsl_matrix* M=gsl_matrix_calloc(n,n);
	while(scanf("%lf", &element)!=EOF && j!=n){
		gsl_matrix_set(M, i, j, element);
		i++;		
		if (i==n){i=0;j++;}
	}
	matrix_print("Matrix M:", M);


	gsl_vector* b=gsl_vector_calloc(n);
	while(scanf("%lf", &element)!=EOF && i!=n){
		gsl_vector_set(b, i, element);
		i++;
	}
	vector_print("Vector b:", b);




	gsl_vector* x=gsl_vector_calloc(n); /* for storing the solution */
	gsl_matrix* Temp=gsl_matrix_calloc(n,n); /* for use in the calculations*/
	gsl_matrix_memcpy(Temp, M);
	gsl_linalg_HH_solve(Temp, b, x);

	vector_print("Solution for x:", x);

	gsl_blas_dgemv(CblasNoTrans, 1, M, x, 0, b);
	printf("Testing solution\n");
	vector_print("calculating Mx (should equal b):", b);

return 0;
}
