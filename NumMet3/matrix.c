#include "matrix.h"
#define RND (int)rand()/RAND_MAX*10

void matrix_generate_random(gsl_matrix* m)
{
	gsl_matrix_set(m, 1, 1, 3);	

	for (size_t j=0; j<m->size2; j++){
		for (size_t i=j; i<m->size1; i++){
			double x = rand()%10;
			gsl_matrix_set(m, i, j, x);
			gsl_matrix_set(m, j, i, x);
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

int matrix_largestElement(gsl_matrix* m, int n)
{	
	assert(n==0 || n==1);
	int index[] = {0,1};
	double largest = gsl_matrix_get(m, 0, 1);
	for (int i=0; i<m->size1; i++){
                for (int j=i+1; j<m->size2; j++){
                        if (gsl_matrix_get(m, i, j) > largest){
				index[0] = i;
				index[1] = j;
				largest = gsl_matrix_get(m, i, j);
			}
                }
        }
	return index[n];
}

