#include <math.h>
#include "lt.h"



#define MYFUNC(x,y,z) ((*func)(x,y,z))


double three_eighths2(double a, double b,double c,double d, int N,double (*func)(double, double, double))
{
    int i = 0;
	double par1 = c;
	double par2 = d;
	double deltax;
	double x = a;
	double efe;
	double integral = 0.0;
	double summ = 0.0;
	
	double peso;
	
	do {    deltax = (b-a)/N;
		
		x = x + deltax;
		
		efe = MYFUNC(x,par1,par2);
		
		
		
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


#undef   MYFUNC(x,y,z) ((*func)(x,y,z))

