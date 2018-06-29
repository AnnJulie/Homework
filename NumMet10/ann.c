#include "ann.h"

ann* ann_alloc(int n, double (*f)(double)){
	ann* network = malloc(sizeof(ann));
	network->n = n;
	network->f = f;
	network->data = gsl_vector_alloc(3*n);
	return network;
}

double ann_feed_forward(ann* network, double x)
{
	double sum = 0;
	for(int i=0; i<network->n; i++){
		double a = gsl_vector_get(network->data, 0*network->n+i);
		double b = gsl_vector_get(network->data, 1*network->n+i);
		double w = gsl_vector_get(network->data, 2*network->n+i);
		sum += network->f((x+a)/b)*w;
	}
	return sum;
}

void ann_train(ann* network, gsl_vector* xlist, gsl_vector* ylist)
{
	double delta(const gsl_vector* p, void* params){
	        gsl_vector_memcpy(network->data, p);
        	double sum = 0;
	        for(int i=0; i<xlist->size; i++){
        	        double x=gsl_vector_get(xlist,i);
                	double f=gsl_vector_get(ylist,i);
	                double y=ann_feed_forward(network, x);
        	        sum+=fabs(y-f);
        	}
        	return sum/xlist->size;
	}

	gsl_vector* p = gsl_vector_alloc(network->data->size);
	gsl_vector_memcpy(p, network->data);

	gsl_vector* step_size = gsl_vector_alloc(network->data->size);
	gsl_vector_set_all(step_size, 0.1);
		
	gsl_multimin_function F;
	F.n = p->size;
	F.f = delta;
	F.params = NULL;

	gsl_multimin_fminimizer *s = 
	gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2, F.n);
	gsl_multimin_fminimizer_set(s, &F, p, step_size);

	int iter = 0;
	int status;
	do{
		iter++;
		int iteration_status = gsl_multimin_fminimizer_iterate(s);
		if(iteration_status != 0){
			printf("System not improving\n");
			break;
		}
		double acc = 0.001;
		status = gsl_multimin_test_size(s->size, acc);
		if(status == GSL_SUCCESS){
			printf("System converged! iterations = %d\n", iter);
		}
	}while(status == GSL_CONTINUE && iter<1000000);
	printf("Network parameters:\n");
	printf("a \t\t b \t\t w \n");
	for (int i=0; i<network->n; i++){
		printf("%g \t %g \t %g \n",
		gsl_vector_get(network->data, 0*network->n + i),
		gsl_vector_get(network->data, 1*network->n + i),
		gsl_vector_get(network->data, 2*network->n + i));	}
	gsl_vector_memcpy(network->data, s->x);
	gsl_vector_free(p);
	gsl_vector_free(step_size);
	gsl_multimin_fminimizer_free(s);
}

void ann_free(ann* network)
{
	gsl_vector_free(network->data);
	free(network);
}




