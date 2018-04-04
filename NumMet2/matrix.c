#include "matrix.h"
#define RND (int)rand()/RAND_MAX*10

void matrix_generate_random(gsl_matrix* m)
{
	gsl_matrix_set(m, 1, 1, 3);	

	for (size_t j=0; j<m->size2; j++){
		for (size_t i=0; i<m->size1; i++){
			gsl_matrix_set(m, i, j, rand()%10);
		}
	}
}

void matrix_print(gsl_matrix* m, char* name)
{	
	printf("Matrix %s:\n", name);
	for (int i=0; i<m->size1; i++){
		for (int j=0; j<m->size2; j++){
			printf("%g ", gsl_matrix_get(m, i, j));
		}
		printf("\n");
	}
	printf("\n");
}



