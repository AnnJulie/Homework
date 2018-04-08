#include "jacobiDiagonalization.h"


// Jacobi diagonalization: upper triangle of A is destroyed;
// e and V accumulate eigenvalues and eigenvectors



int jacobiValueByValue(gsl_matrix* A, gsl_vector* e, gsl_matrix* V)
{
	int changed, sweeps=0, n=A->size1;
        gsl_matrix_set_identity(V);
        for (int i=0; i<n; i++)gsl_vector_set(e, i, gsl_matrix_get(A, i, i));
        do{ changed=0; sweeps++; int p,q;
                p = matrix_largestElement(A,0);
		q = matrix_largestElement(A,1);
	                double app=gsl_vector_get(e,p);
                        double aqq=gsl_vector_get(e,q);
                        double apq=gsl_matrix_get(A,p,q);
                        double phi=0.5*atan2(2*apq, aqq-app);
                        double c = cos(phi), s = sin(phi);
                        double app1=c*c*app - 2*s*c*apq + s*s*aqq;
                        double aqq1=s*s*app + 2*s*c*apq + c*c*aqq;
                        if(app1!=app || aqq1!=aqq){
                                changed = 1;
                                gsl_vector_set(e,p,app1);
                                gsl_vector_set(e,q,aqq1);
                                gsl_matrix_set(A,p,q,0.0);
				for(int i=0;i<p;i++){
                                        double aip = gsl_matrix_get(A,i,p);
                                        double aiq = gsl_matrix_get(A,i,q);
                                        gsl_matrix_set(A,i,p,c*aip-s*aiq);
                                        gsl_matrix_set(A,i,q,c*aiq+s*aip);}
                                for(int i=p+1;i<q;i++){
                                        double api = gsl_matrix_get(A,p,i);
                                        double aiq = gsl_matrix_get(A,i,q);
                                        gsl_matrix_set(A,p,i,c*api-s*aiq);
                                        gsl_matrix_set(A,i,q,c*aiq+s*api);}
                                for(int i=q+1;i<n;i++){
                                        double api = gsl_matrix_get(A,p,i);
                                        double aqi = gsl_matrix_get(A,q,i);
                                        gsl_matrix_set(A,p,i,c*api-s*aqi);
                                        gsl_matrix_set(A,q,i,c*aqi+s*api);}
                                for(int i=0;i<n;i++){
                                        double vip = gsl_matrix_get(V,i,p);
                                        double viq = gsl_matrix_get(V,i,q);
                                        gsl_matrix_set(V,i,p,c*vip-s*viq);
                                        gsl_matrix_set(V,i,q,c*viq+s*vip);}
                                }}while (changed !=0);

        return sweeps;


	
}








int jacobiCyclic(gsl_matrix* A, gsl_vector* e, gsl_matrix* V)
{

	int changed, sweeps=0, n=A->size1;
	gsl_matrix_set_identity(V);
	for (int i=0; i<n; i++)gsl_vector_set(e, i, gsl_matrix_get(A, i, i));
	do{ changed=0; sweeps++; int p,q;
		for(p=0;p<n;p++)for(q=p+1;q<n;q++){
			double app=gsl_vector_get(e,p);
			double aqq=gsl_vector_get(e,q);
			double apq=gsl_matrix_get(A,p,q);
			double phi=0.5*atan2(2*apq, aqq-app);
			double c = cos(phi), s = sin(phi);
			double app1=c*c*app - 2*s*c*apq + s*s*aqq;
			double aqq1=s*s*app + 2*s*c*apq + c*c*aqq;
			if(app1!=app || aqq1!=aqq){ 
				changed = 1;
				gsl_vector_set(e,p,app1);
				gsl_vector_set(e,q,aqq1);
				gsl_matrix_set(A,p,q,0.0);
				for(int i=0;i<p;i++){
					double aip = gsl_matrix_get(A,i,p);
					double aiq = gsl_matrix_get(A,i,q);
					gsl_matrix_set(A,i,p,c*aip-s*aiq);
					gsl_matrix_set(A,i,q,c*aiq+s*aip);}
				for(int i=p+1;i<q;i++){
					double api = gsl_matrix_get(A,p,i);
                                        double aiq = gsl_matrix_get(A,i,q);
                                        gsl_matrix_set(A,p,i,c*api-s*aiq);
                                        gsl_matrix_set(A,i,q,c*aiq+s*api);}
				for(int i=q+1;i<n;i++){
                                        double api = gsl_matrix_get(A,p,i);
                                        double aqi = gsl_matrix_get(A,q,i);
                                        gsl_matrix_set(A,p,i,c*api-s*aqi);
                                        gsl_matrix_set(A,q,i,c*aqi+s*api);}
				for(int i=0;i<n;i++){
                                        double vip = gsl_matrix_get(V,i,p);
                                        double viq = gsl_matrix_get(V,i,q);
                                        gsl_matrix_set(V,i,p,c*vip-s*viq);
                                        gsl_matrix_set(V,i,q,c*viq+s*vip);}
				}}}while (changed !=0);

	return sweeps;
}


int lowestEigenvaluesCyclic(gsl_matrix* A, gsl_vector* e, gsl_matrix* V)
{
	int changed, sweeps=0, eigenvalues=0, n=A->size1;
        gsl_matrix_set_identity(V);
        for (int i=0; i<n; i++)gsl_vector_set(e, i, gsl_matrix_get(A, i, i));
        do{ changed=0; sweeps++; int p,q;
                for(p=0;p<n;p++){
		   for(q=p+1;q<n;q++){
                        double app=gsl_vector_get(e,p);
                        double aqq=gsl_vector_get(e,q);
                        double apq=gsl_matrix_get(A,p,q);
                        double phi=0.5*atan2(2*apq, aqq-app);
                        double c = cos(phi), s = sin(phi);
                        double app1=c*c*app - 2*s*c*apq + s*s*aqq;
                        double aqq1=s*s*app + 2*s*c*apq + c*c*aqq;
                        if(app1!=app || aqq1!=aqq){
                                changed = 1;
                                gsl_vector_set(e,p,app1);
                                gsl_vector_set(e,q,aqq1);
                                gsl_matrix_set(A,p,q,0.0);
                                for(int i=0;i<p;i++){
                                        double aip = gsl_matrix_get(A,i,p);
                                        double aiq = gsl_matrix_get(A,i,q);
                                        gsl_matrix_set(A,i,p,c*aip-s*aiq);
                                        gsl_matrix_set(A,i,q,c*aiq+s*aip);}
                                for(int i=p+1;i<q;i++){
                                        double api = gsl_matrix_get(A,p,i);
                                        double aiq = gsl_matrix_get(A,i,q);
                                        gsl_matrix_set(A,p,i,c*api-s*aiq);
                                        gsl_matrix_set(A,i,q,c*aiq+s*api);}
                                for(int i=q+1;i<n;i++){
                                        double api = gsl_matrix_get(A,p,i);
                                        double aqi = gsl_matrix_get(A,q,i);
                                        gsl_matrix_set(A,p,i,c*api-s*aqi);
                                        gsl_matrix_set(A,q,i,c*aqi+s*api);}
                                for(int i=0;i<n;i++){
                                        double vip = gsl_matrix_get(V,i,p);
                                        double viq = gsl_matrix_get(V,i,q);
                                        gsl_matrix_set(V,i,p,c*vip-s*viq);
                                        gsl_matrix_set(V,i,q,c*viq+s*vip);}
			}
		   }
		   if(p>=1 && gsl_vector_get(e,p) != gsl_vector_get(e, p-1)){
			printf("e[%d] = %g, e[%d] = %g, eigenvalues = %d\n", p-1, gsl_vector_get(e,p-1), p, gsl_vector_get(e,p), eigenvalues);
			return eigenvalues;
		   }else{eigenvalues++;
			printf("p=%d\n",p); }
		}}while (changed !=0);
        return eigenvalues;
}



// generates a random nxn matrix. return diagonalization time
double cyclicDiagTime(size_t n)
{
        clock_t begin = clock();
        // generate random matrix M, allocate vector e and matrix V
        gsl_matrix* M = gsl_matrix_calloc(n,n);
        matrix_generate_random(M);
        gsl_vector* e = gsl_vector_alloc(n);
        gsl_matrix* V = gsl_matrix_alloc(n,n);
        // jacobi diagonalization
        int sweeps = jacobiCyclic(M, e, V);
        clock_t end = clock();
        // free allocated space
        gsl_vector_free(e);
        gsl_matrix_free(M);
        gsl_matrix_free(V);
        return (double)(end - begin) / CLOCKS_PER_SEC;;
}



void testCyclicDiagonalization()
{
        // Generate a 3x3 matrix A with known eigenvalues 0,-3,-3
        size_t n = 3;
        gsl_matrix* A = gsl_matrix_alloc(n,n);
        gsl_matrix_set_all(A,1.0);
        for (size_t i=0; i<n; i++){
                gsl_matrix_set(A,i,i,-2.0);
        }
        matrix_print(A, "A with eigenvalues 0,-3,-3");

        //allocate vector e and matrix V
        gsl_vector* e = gsl_vector_alloc(n);
        gsl_matrix* V = gsl_matrix_alloc(n,n);

        // jacobi diagonalizationd
        int sweeps = jacobiCyclic(A, e, V);
	printf("Diagonalized A with %d sweeps\n", sweeps);
        matrix_print(V, "V");
        vector_print(e, "e");

        //free allocated space
        gsl_vector_free(e);
        gsl_matrix_free(A);
        gsl_matrix_free(V);
}




// generates a random nxn matrix. return diagonalization time
double valueDiagTime(size_t n)
{       
        clock_t begin = clock();
        // generate random matrix M, allocate vector e and matrix V
        gsl_matrix* M = gsl_matrix_calloc(n,n);
        matrix_generate_random(M);
        gsl_vector* e = gsl_vector_alloc(n);
        gsl_matrix* V = gsl_matrix_alloc(n,n);
        // jacobi diagonalization
        int sweeps = jacobiValueByValue(M, e, V);
        clock_t end = clock();
        // free allocated space
        gsl_vector_free(e);
        gsl_matrix_free(M);
        gsl_matrix_free(V); 
        return (double)(end - begin) / CLOCKS_PER_SEC;;
}
