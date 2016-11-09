/*
*                                                           
* FUNCTION: rdf_ipx                                        
*                                                       
* PURPOSE: reads data from a file and convert it in a double area
*                                                      
* INPUT:  file         = output file                 
*         dim          = dimension of the data      
*         nx,ny         = number of data in each dimension
*         type         = type of data (int,double,double)
*         mode         = mode of storage (bin,txt)     
*         nature       = nature of data (real,complex)
*         comments     = file comments               
*                                                   
* VERSION: 2.0  Mars  1994
*                       
* AUTHOR: JP Kneib
*                     
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "lt.h"

/**************************************************************/
double	**rdf_ipx(	char *file,
					int *nx,int *ny,
					char *type,char *mode,char *nature,char comments[1024],
					double *xmin,double *xmax,double *ymin,double *ymax)
{
/****************  declarations  **************************/
FILE		*fp;
register	int	i,j;
auto		int	dim,ival;
auto		double	dval;
auto		int	*vect_i;
auto		double	*vect_d;
auto		double	**square_f_real;

/*****************  verification du format du fichier  *****************/
if ((fp=fopen(file,"r")) != NULL)
        {
        fscanf(fp,"%d%lf%lf%lf%lf\n",&dim,xmin,xmax,ymin,ymax);
        fscanf(fp,"%d\n",nx);
        fscanf(fp,"%d\n",ny);
        fscanf(fp,"%s %s %s\n",type,mode,nature);
        fgets(comments,1024,fp);
	fscanf(fp,"\n");
        }
else
	{
fprintf(stderr,"\n\nFATAL ERROR reading file %s that doesn't exist\n",file);
	exit(-1);
        }

square_f_real = (double **)alloc_square_int(*nx,*ny);

/************  image 2D  ************************/

	if(strcmp(mode,"bin") == 0)
                {
                if ((strcmp(type,"int") == 0) && (strcmp(nature,"real") == 0))
                        {
			vect_i=(int *)malloc(*ny);
			for (i=0;i<*nx;i++)
			{
                        fread(vect_i,sizeof(int),*ny,fp);
			for (j=0;j<*ny;j++)
				square_f_real[i][j]=(double) vect_i[j];
			}
			free(vect_i);
                        }
                if ((strcmp(type,"double") == 0) && (strcmp(nature,"real") == 0))
                        {
			for (i=0;i<*nx;i++)
                        fread(square_f_real[i],sizeof(double),*ny,fp);
                        }
              if ((strcmp(type,"double") == 0) && (strcmp(nature,"real") == 0))
                        {
	                vect_d = (double *)malloc(*ny);
			for (i=0;i<*nx;i++)
			{
                        fread(vect_d,sizeof(double),*ny,fp);
                        for (j=0;j<*ny;j++)
                                square_f_real[i][j]=(double) vect_d[j];
                        }
			free(vect_d);
                        }
		
                }

	else if(strcmp(mode,"txt") == 0)
                {
            if ((strcmp(type,"int") == 0) && (strcmp(nature,"real") == 0))
                        {
                        for (i=0;i<*nx;i++)
			for (j=0;j<*ny;j++)
				{
                        	fscanf(fp,"%d\n",&ival);
				square_f_real[i][j]=(double) ival;
				}
                        }
           if ((strcmp(type,"double") == 0) && (strcmp(nature,"real") == 0))
                        {
                        for (i=0;i<*nx;i++)
			for (j=0;j<*ny;j++)
                        	fscanf(fp,"%lf\n",&square_f_real[i][j]);
                        }
           if ((strcmp(type,"double") == 0) && (strcmp(nature,"real") == 0))
                        {
                        for (i=0;i<*nx;i++)
			for (j=0;j<*ny;j++)
				{
                        	fscanf(fp,"%lf\n",&dval);
                                square_f_real[i][j]=(double) dval;
                                }
                        }

                }
	else
                {
fprintf(stderr,"\n\nFATAL ERROR file %s have unknown mode %s\n",file,mode);
                exit(-1);
                }

return((double **)square_f_real);
}
