#include <math.h>
#include "lt.h"



#define FUNC(x) ((*func)(x))


double three_eighths(double a, double b, int N,double (*func)(double))
{
    int i = 0;
	double deltax;
	double x = a;
	double efe;
	double integral = 0.0;
	double summ = 0.0;
	
	double peso;
	
	do {    deltax = (b-a)/N;
		
		x = x + deltax;
		
		efe = FUNC(x);
		
		
		
		if(i == 0){ peso = (3.0*deltax)/8.0; }	      
		else if ( i == N ){ peso = (3.0*deltax)/8.0 ; }
		else if (i%3 == 0) { peso = (9.0*deltax)/8.0;}
		else if ( (i+1)%3 == 0) { peso = (9.0*deltax)/8.0;}       
		else {peso = (6.0*deltax)/8.0;}
		
		integral = efe*peso;
		
		summ = summ + integral;
		
		
		i = i + 1;
		
		
	} while (    i <= N  );
	
	
	
	
	return(summ);
	

}


#undef   FUNC(x) ((*func)(x))

