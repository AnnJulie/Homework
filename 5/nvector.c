#include "nvector.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

nvector *nvector_alloc(int n){
	nvector* v = malloc(sizeof(nvector));
	(*v).size = n;
	(*v).data = malloc(n*sizeof(double));
	if( v == NULL ) fprintf(stderr,"error in nvector_alloc\n");
	return v;
}

void nvector_free(nvector* v){ free(v->data); free(v); }

void nvector_set(nvector* v, int i, double value){ (*v).data[i] = value; }

double nvector_get(nvector* v, int i){ return (*v).data[i]; }

double nvector_dot_product(nvector* u, nvector* v){
	if( (*u).size != (*v).size ){ 
		fprintf(stderr, "vectors not equal size\n"); 
		return 0;
	}
	
	double dot_product=0;
	int i;
	for (i=0; i<(*u).size; i++){
		dot_product += ( (*u).data[i] * (*v).data[i] );
	}
	return dot_product;	
}

int double_equal(double a, double b){
	if ( fabs(a-b)/2  == fabs(a-b) ) {return 1;}
	else {return 0;}
}

void nvector_print(nvector* v){
	int i = 0;
	for (i=0; i<(*v).size; i++){
		printf("%.1f ", (*v).data[i]);
	}
	printf("\n");

}
