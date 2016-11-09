/*This code will tabulate gnfw profile parameters for look-up in lenstool*/
/*Last modified: 10/18/04*/
//Last modified: 08/27/07 (EJ)

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
/*#include "nr.h"*/
#include "nrutil.h"
#include "lenstool_tab_func.h"

typedef struct{
	double alpha_now,x_now,kappa,dpl;
} lensdata;

int main(int argc, char *argv[])
{
	//important variables
	char infile[20], outfile[200], one_line[200], *dir, *line_ptr;
	double x_not,a_step,x_max,alphalow,alphahigh,alphastep;
	int alphacount,xcount,countsfile;/*for file i/o*/
	double n_xstepsdouble, n_alphastepsdouble;
	long int n_xstepsint, n_alphastepsint;
	long int stepstot,count;
	lensdata lensprop;
		
	FILE *inp;
	FILE *outf;

	//command line switches*/
	printf("Hello - welcome to lenstool_tab\n");
	// printf("No of arguments = %d\n",argc);
	if(argc != 2)
	{
		printf("The current form of the program is: \n");
		printf("  lenstool_tab [name of input parameter file]\n");
		exit(-1);
	}
	else 
	{
		strcpy(infile, argv[1]);
		inp = fopen(infile, "r");
            // make sure we write to the right directory!
            dir = getenv("LENSTOOL_DIR");
         if( dir == 0 ){
             fprintf(stderr, "ERROR: LENSTOOL_DIR environment variable is not defined\n");
             exit(1);
         }

	      strcpy(outfile,dir);
            strcat(outfile,"/table_src/lenstool.tab");
		printf("The output database will be written to\n  %s\n", outfile);
		// strcpy(outfile, "lenstool.tab");
		outf = fopen(outfile, "wb+");
	}

	//INPUT FILE READ
	while (1)
	{
		line_ptr = fgets(one_line,sizeof(one_line),inp);
		if (line_ptr == NULL) break;
		sscanf(one_line,"%lf %lf %lf %lf %lf %lf", 
			&alphalow, &alphahigh, &alphastep, &x_not, &x_max, &a_step);
	} 

	printf ("Input file %s read in OK:\n",infile);
	printf (" %f %f %f %f %f %f\n", alphalow, alphahigh, alphastep, x_not, x_max, a_step);
	fclose(inp);
	
	//Compute total number of steps for later array manipulation
	n_alphastepsdouble = (alphahigh - alphalow)/alphastep + 1.01;
	n_xstepsdouble = log(x_max/x_not)/log(a_step)+1.01;
	n_alphastepsint = n_alphastepsdouble;
	n_xstepsint = n_xstepsdouble;
	stepstot = n_alphastepsint*n_xstepsint;
	printf("Steps in alpha: %i\n",n_alphastepsint);
	printf("Steps in x: %i\n", n_xstepsint);
	printf("Total steps: %i \n", stepstot);

	//First step is to write parameters used to make  this table into the first line of the binary file
	//in this order: alpha_low, alpha_nsteps, alpha_stepsize*/
	lensprop.alpha_now = alphalow;
	lensprop.x_now = n_alphastepsint;
	lensprop.kappa = alphastep;
	lensprop.dpl = 0.0;
      fwrite(&lensprop,sizeof(lensdata),1,outf);
      
      //in this order: xlow, x_nsteps, x_stepsize
	lensprop.alpha_now = x_not;
	lensprop.x_now = n_xstepsint;
	lensprop.kappa = a_step;
	lensprop.dpl = 0.0;
	fwrite(&lensprop,sizeof(lensdata),1,outf);


	//First step: inner slope loop
	count = 0;	
	lensprop.alpha_now = alphalow;
	for(  alphacount = 0; alphacount < n_alphastepsint; alphacount++ )
	{
		/*now for the loop in x*/
		lensprop.x_now = x_not;
		for( xcount = 0; xcount < n_xstepsint; xcount++ )
		{
			lensprop.kappa = nfwg_kappa(lensprop.x_now,1.0,1.0,lensprop.alpha_now);
			//note that I set r_sc = 1.0, kappas = 1.0*/
			lensprop.dpl = lensprop.x_now*nfwg_dpl(lensprop.x_now,1.0,1.0,lensprop.alpha_now);
			//note that I set r_sc = 1.0, kappas = 1.0*/
			/*note also the x_now out in front...this is to get rid 'r' dependenc in nfwg_dpl*/
			/*so entire scaling is 1.0...this must be corrected for in JP's code*/
			/*lensprop.gamma = */
			/*not putting in lensprop.gamma for now*/
			/*now write lensprop structure to binary file*/
			fwrite(&lensprop,sizeof(lensdata),1,outf);
			printf("Progress : alpha %lf/%lf  x %lf/%lf\r", 
				lensprop.alpha_now,alphahigh,lensprop.x_now,x_max);
			fflush(stdout);

			lensprop.x_now *= a_step;
			count++; 
		}
		lensprop.alpha_now += alphastep;
	}
	
	printf("\n");
	printf("Total steps = %i\n", count);
	printf("Alpha count = %i\n", alphacount);
	printf("xcount = %i\n", xcount);
	printf("Table written in lenstool.tab\n");
	//now store elements into array*/
	fclose(outf);
//	countsfile = read_bin(count);

	return 0;
}
