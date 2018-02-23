#ifndef HAVE_NVECTOR_H /* for multiple includes */


typedef struct {int size; double* data;} nvector;

nvector* nvector_alloc		(int n);	/* allocates memory for size-n vector*/
void nvector_free		(nvector* v);			/* frees memory */
void nvector_set		(nvector* v, int i, double value); /* v_i <- value */
double nvector_get		(nvector* v, int i);		 /*returns v_i */
double nvector_dot_product	(nvector* u, nvector* v);	 /* returns dot-product*/

/* optional functions */
int double_equal		(double a, double b);	/* returns 1 if equal*/
void nvector_print		(nvector* v);		/* prints the vector */				

#define HAVE_NVECTOR_H
#endif	
