/***************************************************************/
/*                                                             */
/* FUNCTION: write_data                                        */
/*                                                             */
/* PURPOSE: writes data in a file			       */
/*                                                             */
/* INPUT:  file		= output file                          */
/*	   data		= pointer to data		       */
/*	   dim		= dimension of the data		       */
/*	   size		= number of data in each dimension     */
/*	   type		= type of data (int,double,double)      */
/*	   mode		= stroring mode of data (txt,bin)      */
/*	   nature	= nature of data (real,complex)        */
/*	   comments	= file comments			       */
/*                                                             */
/* VERSION: 1.0  May  1992                                     */
/*                                                             */
/* AUTHOR: Karim BOUYOUCEF                                     */
/*                                                             */
/***************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "lt.h"

/**************************************************************/
void	wrf_ipx(char *file,long int data,
				int dim,int *size,
				char *type,char *mode,char *nature,char *comments,
				double xmin,double xmax,double ymin,double ymax)
{
/****************  declarations  **************************/

FILE		*fp;
register	int	i,j;
auto		int	nb;
auto		int	**square_i_real,**datai;
auto		double	*vector_f_real,**square_f_real,**dataf;

/*****************  verification du format du fichier  *****************/
fp=fopen(file,"w");
if (fp != NULL)
	{
        fprintf(fp,"%d\t%lf %lf %lf %lf\n",dim,xmin,xmax,ymin,ymax);
        for (i=0;i<dim;i++)
                fprintf(fp,"%d\n",size[i]);
        fprintf(fp,"%s %s %s\n",type,mode,nature);
	nb = strlen(comments);
	if (isspace(comments[nb-1]) != 0)
        	comments[nb-1] = '\0';
        fputs(comments,fp);
	fprintf(fp,"\n");
        }
else
	{
        fprintf(stderr,"\n\nFATAL ERROR writing file %s \n",file);
        exit(-1);
        }

if ((strcmp(type,"int") != 0) && (strcmp(type,"double") != 0) )
	{
fprintf(stderr,"\n\nFATAL ERROR file %s must be of type %s\n",file,type);
	exit(-1);
	}
if ((strcmp(mode,"bin") != 0) && (strcmp(mode,"txt") != 0))
	{
	fprintf(stderr,"\n\nFATAL ERROR file %s must be of mode bin or txt\n",file);
	exit(-1);
	}
if ((strcmp(nature,"real") != 0) && (strcmp(nature,"complex") != 0))
	{
fprintf(stderr,"\n\nFATAL ERROR file %s must be of nature real or complex\n",file);
	exit(-1);
	}

/********************  Ecriture dans le fichier  ***************************/
switch(dim)
	{

case 1:	/****************  signal 1D  ***********************************/
if (strcmp(mode,"bin") == 0)
	{
	if ((strcmp(type,"double") == 0) && (strcmp(nature,"real") == 0))
		{
		vector_f_real =(double *) data;
		fwrite(vector_f_real,sizeof(double),size[0],fp);
		}
	}
else if (strcmp(mode,"txt") == 0)
        {
        if ((strcmp(type,"double") == 0) && (strcmp(nature,"real") == 0))
                {
                vector_f_real =(double *) data;
		for (i=0;i<size[0];i++)
			fprintf(fp,"%lf\n",vector_f_real[i]);
                }
	}
	if (strcmp(mode,"bin") == 0)
		fprintf(fp,"\n");
	fclose(fp);
	break;

case 2:	/************  image 2D  ************************/

if (strcmp(mode,"bin") == 0)
        {
        if ((strcmp(type,"int") == 0) && (strcmp(nature,"real") == 0))
                {
		square_i_real =(int **) alloc_square_int(size[1],size[0]);
		datai =(int **) alloc_square_int(size[0],size[1]);
		datai =(int **) data;
		for(i=0;i<size[1];i++)
		for(j=0;j<size[0];j++)
			square_i_real[i][j]=datai[j][i];

		for (i=0;i<size[1];i++)
                	fwrite(square_i_real[i],sizeof(int),size[0],fp);
                }
        if ((strcmp(type,"double") == 0) && (strcmp(nature,"real") == 0))
                {
		square_f_real =(double **) alloc_square_double(size[1],size[0]);
		dataf =(double **) alloc_square_double(size[0],size[1]);
		dataf =(double **) data;
		for(i=0;i<size[1];i++)
		for(j=0;j<size[0];j++)
			square_f_real[i][j]=dataf[j][i];

		for (i=0;i<size[1];i++)
                        fwrite(square_f_real[i],sizeof(double),size[0],fp);
                }
	}
else if (strcmp(mode,"txt") == 0)
        {
        if ((strcmp(type,"int") == 0) && (strcmp(nature,"real") == 0))
                {
                square_i_real =(int **) data;
                for (i=0;i<size[1];i++)
		for (j=0;j<size[0];j++)
                        fprintf(fp,"%d\n",square_i_real[j][i]);
                }
        if ((strcmp(type,"double") == 0) && (strcmp(nature,"real") == 0))
                {
                square_f_real =(double **) data;
                for (i=0;i<size[1];i++)
		for (j=0;j<size[0];j++)
                        fprintf(fp,"%lf\n",square_f_real[j][i]);
                }
        }
	if (strcmp(mode,"bin") == 0)
                fprintf(fp,"\n");
        fclose(fp);
        break;

default : /***************  dimension > 2  ***************/
 	  fprintf(stderr,"\nFATAL ERROR data must be of dimension < 3\n");
          exit(-1);		
	  break;
	}
}
