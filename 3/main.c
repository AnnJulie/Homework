#include "stdio.h"
#include "limits.h"
#include "float.h"

void myMaxInt()
{
	printf("Testing the 'while-loop'...\n");
	int i = 1;
	while (i+1 > i) {i++;}
	printf("my max int = %i\n", i);

	printf("Testing the 'for-loop'..\n");
	i = 1;
	for (i = 1; i+1 > i; i++){}
	printf("my max int = %i\n", i);

	printf("Testing the 'do-while-loop'...\n");
	i = 1;
	do{i++;}while(i+1 > i);
	printf("my max int = %i\n", i);
}


void myEpsilon()
{
	float xf = 1;
	double xd = 1;
	long double xld = 1L;
	
	// float
	while(1+xf!=1){xf/=2;} xf*=2;
	printf("my 'float' epsilon = %f \n", xf);

	//double
	while(1+xd!=1){xd/=2;} xd*=2;
	printf("my 'double' epsilon = %lg \n", xd);

	// long double
	while(1+xld!=1){xld/=2;} xld*=2;
	printf("my 'long double' epsilon = %Lg \n", xld);

}

void mySum()
{

	// When doing this with 'float', the upward and downward sum did not converge to the same sum. This is due to many small rounding errors of the small storage capacity of the float type.
	int max = INT_MAX/4;
	int i = 1;
	double sum_up = 0.0;
	double sum_down = 0.0;

	while (i<=max){
		sum_up += (1.0/i);
		i++;
	}
	printf("The upward sum is: %lg\n", sum_up);

	i=max;
	while (i>0){
		sum_down += (1.0L/i);
		i--;
	}
	printf("The downward sum is: %lg\n", sum_down);
}

int equal(double a, double b, double tau, double epsilon);


int main(){
	
//	myMaxInt();
//	myEpsilon();
//	mySum();
	equal(3.0, 2.0, 1.0, 1.0);
	return 0;
}
