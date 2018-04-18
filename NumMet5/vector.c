#include "vector.h"

// assign a matrix column to a vector
void vector_get_from_matrix(gsl_vector* v, gsl_matrix* m, size_t col)
{
	assert(v->size == m->size1);
	for (int i=0; i<v->size; i++){
		gsl_vector_set(v, i, gsl_matrix_get(m, i, col));
	}
}

// assign a vector to a matrix column
void vector_set_to_matrix(gsl_vector* v, gsl_matrix*m, size_t col)
{
	assert(v->size == m->size1);
	for (int i=0; i<v->size; i++){
		gsl_matrix_set(m, i, col, gsl_vector_get(v, i));
	}
}

// assign the normalized vector v1 to vector v2
void vector_store_normalized(gsl_vector* v1, gsl_vector* v2){
	assert(v1->size == v2->size);
	for (int i=0; i<v1->size; i++){
		gsl_vector_set(v2, i, gsl_vector_get(v1, i)/sqrt(vector_inner_product(v1, v1)));
	}
}

// update v1 as part of the gram-schmidt orthogonalization. aj = aj - <qi|aj> qi,
// where v1 = aj and v2 = qi
void vector_update_orthogonalized(gsl_vector* v1, gsl_vector* v2)
{
	assert(v1->size == v2->size);
	double ip = vector_inner_product(v1, v2);
	for (int i=0; i<v1->size; i++){
		gsl_vector_set(v1, i,gsl_vector_get(v1, i) -  gsl_vector_get(v2, i)*ip);
	}
}

// print vector
void vector_print(gsl_vector* v, char* name){
	printf("Vector %s:\n", name);
        for (int i=0; i<v->size; i++){
		printf("%g ", gsl_vector_get(v, i));
        }
        printf("\n\n");
}

// return <v1|v2>
double vector_inner_product(gsl_vector* v1, gsl_vector* v2)
{
	assert(v1->size == v2->size);
	double res = 0;
	for(int i=0; i<v1->size; i++){
		res += v1->data[i] * v2->data[i];
	}
	return res;
}
