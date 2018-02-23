#include "nvector.h"
#include "nvector.h"
#include <stdio.h>
#include <stdlib.h>
#define RND (double)rand()/RAND_MAX

int main()
{
	int n = 5;

	printf("\nmain: testing nvector_alloc ...\n");
	nvector *v = nvector_alloc(n);
	if (v == NULL) printf("test failed\n");
	else printf("test passed\n");

	printf("\nmain: testing nvector_set and nvector_get ...\n");
	double value = RND;
	int i = n/2;
	nvector_set(v, i, value);
	double vi = nvector_get(v, i);
	if (double_equal(vi, value)) {printf("test passed\n");}
 	else {printf("test failed\n");}


	/* making some vectors */
	printf("Main: testing nvector_print\n");
	nvector *a = nvector_alloc(3);
	nvector_set(a, 0, 1);
	nvector_set(a, 1, 1);
	nvector_set(a, 2, 1);
	nvector *b = nvector_alloc(3);
	nvector_set(b, 0, 2);
        nvector_set(b, 1, 2);
        nvector_set(b, 2, 2);
	nvector *c = nvector_alloc(3);

	nvector_print(a);
	nvector_print(b);

	printf("\nmain: testing nvector_dot_product ...\n");
	double dot_product = nvector_dot_product(a, b);
	printf("a*b = %f\n", dot_product);

	printf("\n");
	return 0;
}
