/*this function will read a binary file containing gnfw slope (alpha),
  x (r/rsc), kappa and dpl.  These values will have to be rescaled
  before being used by lenstool.  The function will return the number
  of complete entries made*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#define ARRAY_SIZE 3000000

typedef struct{
	double alpha_now,x_now,kappa,dpl;
} lensdata;

 
//extern lensdata *lens_table;

int read_bin(int tot_count)
{
	int i=0;
	char outfile[20];
	double alphasteps, xsteps;
	int int_alphasteps, int_xsteps;
	
	lensdata *lens_table;  // comment here if extern lens_table uncommented
	lensdata array_tmp[2];
	FILE *outf;
	
	
	strcpy(outfile,"lenstool.tab");
	outf = fopen(outfile,"rb+");

	//read in the first two lines to find out how many more lines to read*/
	fread(&array_tmp[0], sizeof(lensdata),1,outf);
	alphasteps = array_tmp[0].x_now+0.01;
	int_alphasteps = alphasteps;
	
	fread(&array_tmp[1], sizeof(lensdata),1,outf); 
	xsteps = array_tmp[1].x_now+0.01;
	int_xsteps = xsteps;
	
	tot_count = int_xsteps*int_alphasteps + 2;
	
	// Allocate lens_table
	lens_table = (lensdata*) malloc( (unsigned int) tot_count*sizeof(lensdata) );
	
	// Fill the first 2 rows of lens_table with array_tmp
	memcpy( &lens_table[0], &array_tmp[0], sizeof(lensdata) );
	memcpy( &lens_table[1], &array_tmp[1], sizeof(lensdata) );

	if(tot_count > ARRAY_SIZE)
	{
		printf("Need to increase output array size!!!\n");
		return -1;
	}

	for (i=2; i<tot_count; i++)
	{
		fread(&lens_table[i], sizeof(lensdata),1,outf);
		/*printf("lens_out[%i].alpha=%f\n",i,lens_out[i].alpha_now);*/
	}
	printf("Finished reading binary file into array.\n");
	fclose(outf);

	return tot_count;
}
