/***************************************************************/
/*                                                             */
/* FUNCTION: read_data                                         */
/*                                                             */
/* PURPOSE: reads data from a file			       */
/*                                                             */
/* INPUT:  file         = output file                          */
/*         dim          = dimension of the data                */
/*         size         = number of data in each dimension     */
/*         type         = type of data (int,double,double)      */
/*         mode         = mode of storage (bin,txt)            */
/*         nature       = nature of data (real,complex)        */
/*         comments     = file comments                        */
/*                                                             */
/* VERSION: 1.0  May  1992                                     */
/*                                                             */
/* AUTHOR: Karim BOUYOUCEF                                     */
/*                                                             */
/***************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "lt.h"

/**************************************************************/
long int	rd_ipxs(char *file,
					int dim_r,
					int size[4],
					char *type_r,
					char mode[4],
					char *nature_r,
					char comments[1024])
{
/****************  declarations  **************************/
FILE		*fp;
register	int	i,j;
auto		int	dim;
auto		char	type[8];
auto		char	nature[8];
auto		int		**square_i_real;
auto		double	*vector_f_real,**square_f_real;

/*****************  verification du format du fichier  *****************/
if ((fp=fopen(file,"r")) != NULL)
        {
        fscanf(fp,"%d\n",&dim);
        for (i=0;i<dim;i++)
                fscanf(fp,"%d\n",&size[i]);
        fscanf(fp,"%s %s %s\n",type,mode,nature);
        fgets(comments,1024,fp);
	fscanf(fp,"\n");
        }
else
	{
fprintf(stderr,"\n\nFATAL ERROR reading file %s that doesn't exist\n",file);
	exit(-1);
        }
fclose(fp);

if (dim_r != dim)
	{
fprintf(stderr,"\nFATAL ERROR file %s must be of dimension %d\n",file,dim_r);
	exit(-1);
	}
if (strcmp(type,type_r) != 0)
	{
fprintf(stderr,"\n\nFATAL ERROR file %s must be of type %s\n",file,type_r);
	exit(-1);
	}
if (strcmp(nature,nature_r) != 0)
	{
fprintf(stderr,"\n\nFATAL ERROR file %s must be of natue %s\n",file,nature_r);
	exit(-1);
	}

/********************  lecture du fichier  ***************************/
fp=fopen(file,"r");
fscanf(fp,"%d\n",&dim);
for (i=0;i<dim;i++)
        fscanf(fp,"%d\n",&size[i]);
fscanf(fp,"%s %s %s\n",type,mode,nature);
fgets(comments,1024,fp);

switch(dim)
	{

case 1:	/****************  signal 1D  ***********************************/
	if(strcmp(mode,"bin") == 0)
		{
		if ((strcmp(type,"double") == 0) && (strcmp(nature,"real") ==0))
			{
			vector_f_real = (double *)malloc(size[0]*sizeof(double));
			fread(vector_f_real,sizeof(double),size[0],fp);
			return((long int)vector_f_real);
			}
		}

	else if (strcmp(mode,"txt") == 0)
		{
                if ((strcmp(type,"double") == 0) && (strcmp(nature,"real") == 0))
                        {
                        vector_f_real = (double *)malloc(size[0]*sizeof(double));
			for (i=0;i<size[0];i++)
                        	fscanf(fp,"%lf\n",&vector_f_real[i]);
			return((long int)vector_f_real);
                        }
		}	
	else 
		{
fprintf(stderr,"\n\nFATAL ERROR file %s have unknown mode %s\n",file,mode);
        	exit(-1);
		}	
		break;

case 2:	/************  image 2D  ************************/
	if(strcmp(mode,"bin") == 0)
                {
                if ((strcmp(type,"int") == 0) && (strcmp(nature,"real") == 0))
                        {
                        square_i_real = (int **)alloc_square_int(size[0],size[1]);
			for (i=0;i<size[0];i++)
                        fread(square_i_real[i],sizeof(int),size[1],fp);
                        return((long int)square_i_real);
                        }
                if ((strcmp(type,"double") == 0) && (strcmp(nature,"real") == 0))
                        {
                        square_f_real = (double **)alloc_square_double(size[0],size[1]);
			for (i=0;i<size[0];i++)
                        fread(square_f_real[i],sizeof(double),size[1],fp);
                        return((long int)square_f_real);
                        }
                }

	else if(strcmp(mode,"txt") == 0)
                {
            if ((strcmp(type,"int") == 0) && (strcmp(nature,"real") == 0))
                        {
                        square_i_real = (int **)alloc_square_int(size[0],size[1]);
                        for (i=0;i<size[0];i++)
			for (j=0;j<size[1];j++)
                        	fscanf(fp,"%d\n",&square_i_real[i][j]);
                        return((long int)square_i_real);
                        }
           if ((strcmp(type,"double") == 0) && (strcmp(nature,"real") == 0))
                        {
                        square_f_real = (double **)alloc_square_double(size[0],size[1]);
                        for (i=0;i<size[0];i++)
			for (j=0;j<size[1];j++)
                        	fscanf(fp,"%lf\n",&square_f_real[i][j]);
                        return((long int)square_f_real);
                        }
                }
	else
    {
		fprintf(stderr,"\n\nFATAL ERROR file %s have unknown mode %s\n",file,mode);
		exit(-1);
	}
	break;

	default : /***************  dimension > 2  ***************/
		fprintf(stderr,"\nFATAL ERROR file %s must be of dimension %d\n",file,dim_r);
        exit(-1);		
		break;

	}
	return 0;
}
